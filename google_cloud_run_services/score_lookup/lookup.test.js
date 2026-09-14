/*
 * Tests for the score-lookup function. Run with `npm test`.
 *
 * The "reference data tables" tests query the real tables in gs://spliceai-lookup-reference-data over
 * the network. Set SKIP_NETWORK_TESTS=1 to run only the offline tests.
 */
import assert from 'node:assert/strict'
import { describe, test } from 'node:test'

import { handleRequest } from './index.js'
import { createTables, lookupScores, parseVariantParams } from './lookup.js'

describe('parseVariantParams', () => {
    test('normalizes a valid variant', () => {
        assert.deepEqual(
            parseVariantParams({ hg: '38', chrom: '8', pos: '140300616', ref: 't', alt: 'g' }),
            { variant: { hg: '38', chrom: 'chr8', pos: 140300616, ref: 'T', alt: 'G' } })
        assert.equal(parseVariantParams({ hg: '37', chrom: 'chrX', pos: '1', ref: 'A', alt: 'C' }).variant.chrom, 'chrX')
        assert.equal(parseVariantParams({ hg: '38', chrom: 'MT', pos: '100', ref: 'A', alt: 'C' }).variant.chrom, 'chrM')
    })

    test('rejects invalid parameters', () => {
        const valid = { hg: '38', chrom: '8', pos: '140300616', ref: 'T', alt: 'G' }
        for (const [name, value] of [
            ['hg', '19'], ['hg', undefined], ['chrom', '23'], ['chrom', 'chr1_KI270706v1_random'], ['chrom', ''],
            ['pos', '0'], ['pos', '-5'], ['pos', '12a'], ['pos', '1.5'], ['ref', 'X'], ['ref', ''], ['alt', 'A<'],
            ['pos', ['140300616']], ['ref', ['T']],
        ]) {
            const result = parseVariantParams({ ...valid, [name]: value })
            assert.ok(result.error, `expected an error for ${name}=${JSON.stringify(value)}`)
            assert.equal(result.variant, undefined)
        }
    })
})

/* A stand-in for TabixIndexedFile that returns the given rows, or throws the given error. */
const fakeFile = (rows, error) => ({
    getLines: async (chrom, start, end, { lineCallback }) => {
        if (error) {
            throw error
        }
        rows.forEach((row) => lineCallback(row.join('\t')))
    },
})

/* The real table definitions for one genome version, with each table's file swapped for a fake. */
const fakeTables = (hg, rowsOrErrorByTableName) => ({
    [hg]: createTables()[hg].map((table) => {
        const rowsOrError = rowsOrErrorByTableName[table.name] || []
        return { ...table, file: rowsOrError instanceof Error ? fakeFile([], rowsOrError) : fakeFile(rowsOrError) }
    }),
})

const VARIANT = { hg: '38', chrom: 'chr8', pos: 100, ref: 'T', alt: 'G' }

describe('lookupScores', () => {
    test('returns the scores from rows that match the variant exactly', async () => {
        const tables = fakeTables('38', {
            'PrimateAI-3D and PromoterAI': [
                ['chr8', '100', 'T', 'C', '0.9', '0.8', '0.7'],  // different ALT
                ['chr8', '101', 'T', 'G', '0.9', '0.8', '0.7'],  // different position
                ['chr8', '100', 'T', 'G', '0.343', '0.71', '-0.084'],
            ],
            'AlphaMissense': [['chr8', '100', 'T', 'G', 'hg38', 'Q8NH21', 'ENST00000335137.4', 'V2L', '0.2937', 'likely_benign']],
            'AlphaGenome AVI': [['chr8', '100', 'T', 'G', '1.590', '30.48']],
        })
        assert.deepEqual(await lookupScores(tables, VARIANT), {
            scores: {
                primateai3d: { percentile: 0.343, genePercentileThreshold: 0.71 },
                promoterai: { score: -0.084 },
                alphamissense: { score: 0.2937 },
                alphagenome_avi: { rawScore: 1.59, phredScore: 30.48 },
            },
            errors: [],
        })
    })

    test('leaves out scores whose values are missing', async () => {
        const tables = fakeTables('38', {
            'PrimateAI-3D and PromoterAI': [['chr8', '100', 'T', 'G', '', '', '-0.002']],
            'AlphaGenome AVI': [['chr8', '100', 'T', 'G', '', '']],
        })
        assert.deepEqual(await lookupScores(tables, VARIANT), { scores: { promoterai: { score: -0.002 } }, errors: [] })
    })

    test('reports a failed table without losing the other tables\' scores', async () => {
        const tables = fakeTables('38', {
            'AlphaMissense': new Error('HTTP 503'),
            'AlphaGenome AVI': [['chr8', '100', 'T', 'G', '1.590', '30.48']],
        })
        assert.deepEqual(await lookupScores(tables, VARIANT), {
            scores: { alphagenome_avi: { rawScore: 1.59, phredScore: 30.48 } },
            errors: [{ table: 'AlphaMissense', message: 'HTTP 503' }],
        })
    })
})

/* A minimal Express-style response that records what the handler sent. */
const fakeResponse = () => {
    const res = { headers: {}, statusCode: undefined, body: undefined }
    res.set = (name, value) => { res.headers[name] = value; return res }
    res.status = (code) => { res.statusCode = code; return res }
    res.json = (body) => { res.body = body; return res }
    res.send = (body) => { res.body = body; return res }
    return res
}

describe('handleRequest', () => {
    const query = { hg: '38', chrom: '8', pos: '100', ref: 'T', alt: 'G' }

    test('answers a valid GET with the scores and a cacheable, cross-origin response', async () => {
        const res = fakeResponse()
        const tables = fakeTables('38', { 'AlphaGenome AVI': [['chr8', '100', 'T', 'G', '1.590', '30.48']] })
        await handleRequest({ method: 'GET', query }, res, tables)
        assert.equal(res.statusCode, 200)
        assert.equal(res.headers['Access-Control-Allow-Origin'], '*')
        assert.match(res.headers['Cache-Control'], /max-age=\d+/)
        assert.deepEqual(res.body, {
            variant: { hg: '38', chrom: 'chr8', pos: 100, ref: 'T', alt: 'G' },
            scores: { alphagenome_avi: { rawScore: 1.59, phredScore: 30.48 } },
            errors: [],
        })
    })

    test('does not let a response with a failed table be cached', async () => {
        const res = fakeResponse()
        await handleRequest({ method: 'GET', query }, res, fakeTables('38', { 'AlphaMissense': new Error('timeout') }))
        assert.equal(res.statusCode, 200)
        assert.equal(res.headers['Cache-Control'], 'no-store')
        assert.equal(res.body.errors.length, 1)
    })

    test('rejects invalid parameters with a 400', async () => {
        const res = fakeResponse()
        await handleRequest({ method: 'GET', query: { ...query, pos: 'abc' } }, res, fakeTables('38', {}))
        assert.equal(res.statusCode, 400)
        assert.ok(res.body.error)
    })

    test('answers a CORS preflight and rejects other methods', async () => {
        const preflight = fakeResponse()
        await handleRequest({ method: 'OPTIONS', query: {} }, preflight, fakeTables('38', {}))
        assert.equal(preflight.statusCode, 204)
        assert.equal(preflight.headers['Access-Control-Allow-Origin'], '*')

        const post = fakeResponse()
        await handleRequest({ method: 'POST', query }, post, fakeTables('38', {}))
        assert.equal(post.statusCode, 405)
    })
})

describe('reference data tables', { skip: process.env.SKIP_NETWORK_TESTS === '1' }, () => {
    const tables = createTables()
    const lookup = (hg, chrom, pos, ref, alt) => lookupScores(tables, parseVariantParams({ hg, chrom, pos, ref, alt }).variant)

    test('hg38 PrimateAI-3D and PromoterAI', async () => {
        const { scores, errors } = await lookup('38', '7', '117480098', 'C', 'A')
        assert.deepEqual(errors, [])
        assert.deepEqual(scores.primateai3d, { percentile: 0.343, genePercentileThreshold: 0.71 })
        assert.deepEqual(scores.promoterai, { score: 0.084 })
    })

    test('hg19 PrimateAI-3D and PromoterAI', async () => {
        const { scores, errors } = await lookup('37', '7', '117513390', 'G', 'A')
        assert.deepEqual(errors, [])
        assert.deepEqual(scores.primateai3d, { percentile: 0.238, genePercentileThreshold: 0.82 })
        assert.deepEqual(scores.promoterai, { score: -0.023 })
    })

    test('AlphaMissense on both genome versions', async () => {
        for (const hg of ['37', '38']) {
            const { scores, errors } = await lookup(hg, 'chr1', '69094', 'G', 'T')
            assert.deepEqual(errors, [])
            assert.deepEqual(scores.alphamissense, { score: 0.2937 }, `hg${hg}`)
        }
    })

    test('AlphaGenome AVI on hg38 only', async () => {
        const hg38 = await lookup('38', '8', '140300616', 'T', 'G')
        assert.deepEqual(hg38.errors, [])
        assert.deepEqual(hg38.scores.alphagenome_avi, { rawScore: 1.59, phredScore: 30.48 })

        const hg19 = await lookup('37', '8', '140300616', 'T', 'G')
        assert.deepEqual(hg19.errors, [])
        assert.equal(hg19.scores.alphagenome_avi, undefined)
    })

    test('an insertion matches no rows', async () => {
        const { scores, errors } = await lookup('38', '1', '1042601', 'A', 'AGAGAG')
        assert.deepEqual({ scores, errors }, { scores: {}, errors: [] })
    })
})
