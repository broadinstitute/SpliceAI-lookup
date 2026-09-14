/*
 * Looks up one variant's AlphaMissense, PrimateAI-3D, PromoterAI, and AlphaGenome Variant Impact (AVI)
 * scores in the tabix-indexed tables in gs://spliceai-lookup-reference-data.
 *
 * The page used to query these tables from the browser, which downloads each table's whole .tbi index
 * (about 4 MiB for hg38) on every new session and is billed as internet egress. Running the query here,
 * in the bucket's region, means the only data sent to the browser is a few hundred bytes of JSON.
 */
import { TabixIndexedFile } from '@gmod/tabix'

const REFERENCE_DATA_URL = 'https://storage.googleapis.com/spliceai-lookup-reference-data'

// How long one table query may take before it is reported as failed, so a stalled read of one
// table can't hold up the response for the others.
const TABLE_QUERY_TIMEOUT_MS = 10000

/* Each table's columns start with CHROM, POS, REF, ALT. parseRow receives the fields of a row that
 * matches the variant and returns the scores it holds, keyed the way the page's 'Other scores'
 * table expects, leaving out any score whose value is missing. */
const parsePrimateAiAndPromoterAiRow = (fields) => {
    // chrom, pos, ref, alt, PAI3D_percentile, PAI3D_gene_threshold, PromoterAI_score
    const scores = {}
    const percentile = parseFloat(fields[4])
    const genePercentileThreshold = parseFloat(fields[5])
    if (!isNaN(percentile) && !isNaN(genePercentileThreshold)) {
        scores.primateai3d = { percentile, genePercentileThreshold }
    }
    const promoterAiScore = parseFloat(fields[6])
    if (!isNaN(promoterAiScore)) {
        scores.promoterai = { score: promoterAiScore }
    }
    return scores
}

const parseAlphaMissenseRow = (fields) => {
    // CHROM, POS, REF, ALT, genome, uniprot_id, transcript_id, protein_variant, am_pathogenicity, am_class
    const score = parseFloat(fields[8])
    return isNaN(score) ? {} : { alphamissense: { score } }
}

const parseAlphaGenomeAviRow = (fields) => {
    // #CHROM, POS, REF, ALT, raw_score, PHRED
    const rawScore = parseFloat(fields[4])
    const phredScore = parseFloat(fields[5])
    return isNaN(rawScore) || isNaN(phredScore) ? {} : { alphagenome_avi: { rawScore, phredScore } }
}

// The library's default cache keeps up to 1 GiB of decompressed chunks per table for 3 minutes, sized
// for a genome browser panning over a region. A lookup reads about 1.4 MiB per table and rarely
// repeats, so with the default, memory grows with every distinct lookup until it outgrows the instance.
const CHUNK_CACHE_BYTES = 16 * 2 ** 20

const tableFile = (fileName) => {
    const url = `${REFERENCE_DATA_URL}/${fileName}`
    return new TabixIndexedFile({ url, tbiUrl: `${url}.tbi`, chunkCacheSize: CHUNK_CACHE_BYTES })
}

/* The tables to query for each genome version. Google DeepMind only released AVI scores for hg38. */
export const createTables = () => ({
    '37': [
        { name: 'PrimateAI-3D and PromoterAI', file: tableFile('PrimateAI_and_PromoterAI_scores.hg19.20250627.tsv.gz'), parseRow: parsePrimateAiAndPromoterAiRow },
        { name: 'AlphaMissense', file: tableFile('AlphaMissense_hg19.tsv.gz'), parseRow: parseAlphaMissenseRow },
    ],
    '38': [
        { name: 'PrimateAI-3D and PromoterAI', file: tableFile('PrimateAI_and_PromoterAI_scores.hg38.20250627.tsv.gz'), parseRow: parsePrimateAiAndPromoterAiRow },
        { name: 'AlphaMissense', file: tableFile('AlphaMissense_hg38.tsv.gz'), parseRow: parseAlphaMissenseRow },
        { name: 'AlphaGenome AVI', file: tableFile('AlphaGenome_AVI_SNV_scores.hg38.tsv.gz'), parseRow: parseAlphaGenomeAviRow },
    ],
})

const CHROM_REGEX = /^(?:chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$/i
const POS_REGEX = /^[1-9][0-9]{0,9}$/
const ALLELE_REGEX = /^[ACGTN]{1,1000}$/i

/* Validates the query parameters and returns { variant } with the chromosome in the tables' chr-prefixed
 * form and the alleles upper-cased, or { error } describing the first invalid parameter. */
export const parseVariantParams = (params) => {
    // A repeated or bracketed parameter (?pos=1&pos=2, ?ref[]=A) arrives as an array or object
    // rather than a string, so anything that isn't a string counts as missing.
    const param = (name) => (typeof params[name] === 'string' ? params[name] : '')
    const hg = param('hg')
    if (hg !== '37' && hg !== '38') {
        return { error: "hg must be '37' or '38'" }
    }
    const chromMatch = CHROM_REGEX.exec(param('chrom'))
    if (!chromMatch) {
        return { error: 'chrom must be 1-22, X, Y, or M, with or without a chr prefix' }
    }
    const pos = param('pos')
    if (!POS_REGEX.test(pos)) {
        return { error: 'pos must be a positive integer' }
    }
    const ref = param('ref')
    const alt = param('alt')
    if (!ALLELE_REGEX.test(ref) || !ALLELE_REGEX.test(alt)) {
        return { error: 'ref and alt must be DNA sequences of A, C, G, T, or N' }
    }
    const chromName = chromMatch[1].toUpperCase()
    return {
        variant: {
            hg,
            chrom: `chr${chromName === 'MT' ? 'M' : chromName}`,
            pos: parseInt(pos, 10),
            ref: ref.toUpperCase(),
            alt: alt.toUpperCase(),
        },
    }
}

/* Queries every table for the variant's genome version in parallel. Returns { scores, errors }: the
 * scores found across all tables, and one { table, message } entry for each table whose query failed,
 * so a failed read of one table costs only that table's scores. */
export const lookupScores = async (tables, variant) => {
    const { hg, chrom, pos, ref, alt } = variant
    const scores = {}
    const errors = []
    await Promise.all(tables[hg].map(async (table) => {
        try {
            await table.file.getLines(chrom, pos - 1, pos, {
                signal: AbortSignal.timeout(TABLE_QUERY_TIMEOUT_MS),
                lineCallback: (line) => {
                    const fields = line.split('\t')
                    if (parseInt(fields[1], 10) === pos && fields[2] === ref && fields[3] === alt) {
                        Object.assign(scores, table.parseRow(fields))
                    }
                },
            })
        } catch (e) {
            errors.push({ table: table.name, message: String(e && e.message ? e.message : e) })
        }
    }))
    return { scores, errors }
}
