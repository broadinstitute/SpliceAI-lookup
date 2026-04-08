// store user settings at the time of the last 'submit', as well as the resulting responses from the server
let lastGenomeVersion, lastSpliceaiResponseJson, lastPangolinResponseJson
let batchVariantResults = []
let currentBatchIndex = 0
let lastBatchFormOptions = null
const BATCH_VARIANT_MAX = 50

// Define constants 
const GENCODE_VERSION = "v49"

const VARIANT_RE = new RegExp(
    "^[\\s]*" +
    "(chr)?([0-9XYMTt]{1,2})" +
    "[-\\p{Pd}\\s:]+" +
    "([0-9,]+)" +
    "[-\\p{Pd}\\s:]*" +
    "([ACGT]+)" +
    "[-\\p{Pd}\\s:>]+" +
    "([ACGT]+)",
    'iu'
)

const TRANSCRIPT_PRIORITY = {
    "MS": 3,
    "MP": 2,
    "C": 1,
    "N": 0,
}

const nameMapForPredictorScores = {
    "primateai3d": "PrimateAI-3D",
    "promoterai": "PromoterAI",
    "cadd": "CADD",
    "revel": "REVEL",
    "revel_max": "REVEL",
    "sift_max": "SIFT (max)",
    "phylop": "PhyloP",
    //"gnomad": "gnomAD",
    "polyphen_max": "PolyPhen (max)",
    "alphamissense": "AlphaMissense",
}


const predictorPointsToColor = {
    "+8": "#ff8177",  // Very Strong Pathogenic (red)
    "+4": "#ffd7d1",  // Strong Pathogenic (light red)
    "+3": "#ffdcb8",  // Pathogenic (light orange-red)
    "+2": "#ffebc5",  // Moderate Pathogenic (light orange)
    "+1": "#fff5dc",  // Supporting Pathogenic (light yellow)
    "0": "#f0f0f0",   // Indeterminate (very light gray)
    "-1": "#d8f3e8",  // Supporting Benign (very light surf green)
    "-2": "#c4eddd",  // Moderate Benign (pale surf green)
    "-3": "#aee8d2",  // Benign (light surf green)
    "-4": "#9ae2c6",   // Strong Benign (soft surf green)
    "-8": "#6fd9ab",   // Very Strong Benign (green)
}

const predictorPointsToLabel = {
    "+8":  "Very Strong Pathogenic",
    "+4":  "Strong Pathogenic",
    "+3":  "Pathogenic",
    "+2":  "Moderate Pathogenic",
    "+1":  "Supporting Pathogenic",
    "0":  "Indeterminate",
    "-1": "Supporting Benign",
    "-2": "Moderate Benign",
    "-3": "Benign",
    "-4": "Strong Benign",
    "-8": "Very Strong Benign",
}

const predictorScoreToPoints = {
    "cadd": (record) => {
        const score = parseFloat(record.score)
        if (score >= 28.1) {            //moderate (pathogenic)
            return "+2"
        } else if (score >= 25.3) {     //supporting (pathogenic)
            return "+1"
        } else if (score > 22.7) {      //indeterminate
            return "0"
        } else if (score > 17.3) {      // supporting (benign)
            return "-1"
        } else if (score > 0.15) {      // moderate (benign)
            return "-2"
        } else {                        // strong (benign)
            return "-4"
        }
    },
    "revel_max": (record) => {
        const score = parseFloat(record.score)
        if (score >= 0.932) {           // strong (pathogenic)
            return "+4"
        } else if (score >= 0.773) {    // moderate (pathogenic)
            return "+2"
        } else if (score >= 0.644) {    // supporting (pathogenic)
            return "+1"
        } else if (score > 0.290) {     // indeterminate
            return "0"
        } else if (score > 0.183) {     // supporting (benign)
            return "-1"
        } else if (score > 0.16) {      // moderate (benign)
            return "-2"
        } else if (score > 0.003) {      // strong (benign)
            return "-4"
        } else {                        // very strong (benign)
            return "-8"
        }
    },
    "sift_max": (record) => {
        const score = parseFloat(record.score)
        if (score <= 0) {   // moderate (pathogenic)
            return "+2"
        } else if (score < 0.001) { // supporting (pathogenic)
            return "+1"
        } else if (score < 0.08) { // indeterminate
            return "0"
        } else if (score < 0.327) { // supporting (benign)
            return "-1"
        } else { // moderate (benign)
            return "-2"
        }
    },
    "polyphen_max": (record) => {
        const score = parseFloat(record.score)
        if (score >= 0.999) {  // moderate (pathogenic)
            return "+2"
        } else if (score >= 0.978) {  // supporting (pathogenic)
            return "+1"
        } else if (score > 0.113) { // indeterminate
            return "0"
        } else if (score > 0.009) {  // supporting (benign)
            return "-1"
        } else {  // moderate (benign)
            return "-2"
        }
    },
    "phylop": (record) => {
        const score = parseFloat(record.score)
        if (score >= 9.741) {    //moderate (pathogenic)
            return "+2"
        } else if (score >= 7.367) {  //supporting (pathogenic)
            return "+1"
        } else if (score > 1.879) {  //indeterminate
            return "0"
        } else if (score > 0.021) {  //supporting (benign)
            return "-1"
        } else { // moderate (benign)
            return "-2"
        }
    },
    "alphamissense": (record) => {
        const score = parseFloat(record.score)
        if (score >= 0.99) {     // +4 (strong)
            return "+4"
        } else if (score >= 0.972) {  // +3
            return "+3"
        } else if (score >= 0.906) {  // +2 (moderate)
            return "+2"
        } else if (score >= 0.792) {  // +1 (supporting)
            return "+1"
        } else if (score >= 0.170) {   // indeterminate
            return "0"
        } else if (score >= 0.1) {  // supporting (benign)
            return "-1"
        } else if (score >= 0.071) {  // moderate (benign)
            return "-2"
        } else {  //  -3 (benign)
            return "-3"
        }
    },
}


//seqr thresholds: https://github.com/broadinstitute/seqr/blob/master/ui/shared/utils/constants.js#L1493-L1531
//gnomAD browser thresholds: https://github.com/broadinstitute/gnomad-browser/blob/main/browser/src/VariantPage/VariantInSilicoPredictors.tsx#L22-L60
const colorMapForPredictorScores = {
    "primateai3d": (record) => {
        if (record.percentile >= record.genePercentileThreshold) { // + 0.1) {
            return "#fccfb8"  // light red
        //} else if (record.percentile >= record.genePercentileThreshold - 0.1) {
        //    return "#fff19d" // light yellow
        } else {
            return "#ffffff"  // white
        }
    },
    "promoterai": (record) => {
        const absScore = Math.abs(parseFloat(record.score))
        if (absScore >= 0.5) {
            return "#fccfb8"    // light red
        } else if (absScore >= 0.1) {
            return "#fff19d"    // light yellow
        } else {
            return "#ffffff"   // white
        }
    },
    "cadd": (record) => predictorPointsToColor[predictorScoreToPoints["cadd"](record)],
    "revel_max": (record) => predictorPointsToColor[predictorScoreToPoints["revel_max"](record)],
    "sift_max": (record) => predictorPointsToColor[predictorScoreToPoints["sift_max"](record)],
    "polyphen_max": (record) => predictorPointsToColor[predictorScoreToPoints["polyphen_max"](record)],
    "phylop": (record) => predictorPointsToColor[predictorScoreToPoints["phylop"](record)],
    "alphamissense": (record) => predictorPointsToColor[predictorScoreToPoints["alphamissense"](record)],
}

const computeHelpTextForPredictorScores = (predictor, record) => {
    // example "This CADD score is in the 'moderate pathogenic' range (+2 points) with a score of 25.3"
    const points = predictorScoreToPoints[predictor](record)
    const label = predictorPointsToLabel[points]
    const predictorLabel = nameMapForPredictorScores[predictor] || predictor
    let helpText = `The ${record.score} ${nameMapForPredictorScores[predictor] || predictor} score is in the ${label} range and would count for ${points} points based on thresholds and points established for missense variants in <a href='https://www.biorxiv.org/content/10.1101/2024.09.17.611902v1.full' target='_blank'>Bergquist et al. 2024</a> and <a href='https://pmc.ncbi.nlm.nih.gov/articles/PMC9748256/' target='_blank'>Pejaver et al. 2022</a>`
    // generate a color legend as an html table:
    helpText += `<br /><br /><b>Legend:</b><br /><br />`
    helpText += `<table style='border: 1px !important; border-collapse: collapse;'>`
    for (const points of ["+8", "+4", "+3", "+2", "+1", "0", "-1", "-2", "-3", "-4", "-8"]) {
        helpText += `<tr style='background-color:${predictorPointsToColor[points]}'><td style='padding:5px;margin:5px;border-radius:0px;min-width:170px'>${predictorPointsToLabel[points]}</td><td style='min-width:100px; text-align:right'>${points} points</td></tr>`
    }
    helpText += `</table>`
    return helpText
}

// TODO: Is this used?
const helpTextForPredictorScores = {
    "primateai3d": (record) => `<a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10187174/' target='_blank'>PrimateAI-3D</a> gene-specific thresholds are provided by the authors. Values above the threshold are shown with a red background to indicate a 'likely deleterious' prediction.`,
    "promoterai": (record) => `<a href='https://www.science.org/doi/10.1126/science.ads7373' target='_blank'>PromoterAI</a>  &nbsp; <a href='https://github.com/Illumina/PromoterAI' target='_blank'><i class='github square icon'></i></a> &nbsp; scores range from -1 to 1 with 0 meaning no activity. Negative values represent under-expression and positive values represent over-expression. A threshold of +/-0.1 is used for high sensitivity, and +/-0.5 for high precision.`,
    "cadd": (record) => computeHelpTextForPredictorScores("cadd", record),
    "revel_max": (record) => computeHelpTextForPredictorScores("revel_max", record),
    "sift_max": (record) => computeHelpTextForPredictorScores("sift_max", record),
    "polyphen_max": (record) => computeHelpTextForPredictorScores("polyphen_max", record),
    "phylop": (record) => computeHelpTextForPredictorScores("phylop", record),
    "alphamissense": (record) => computeHelpTextForPredictorScores("alphamissense", record),
}

const formatScore = (score) => {
    if (score == null) {
        return ""
    }
    return Math.abs(parseFloat(score)).toFixed(2)
}

const getScoreStyle = (score) => {
    score = parseFloat(score)

    if (Math.abs(score) < 0.01) {
        return `style='color:#BBBBBB;'`
    }

    let color
    if (Math.abs(score) >= 0.8) {
        color = "#fccfb8"
    } else if (Math.abs(score) >= 0.5) {
        color = "#fff19d"
    } else if (Math.abs(score) >= 0.2) {
        color = "#cdffd7"
    } else {
        return ""
    }

    return `style='white-space:nowrap;background-color:${color};'`
}

// Preparing TabixIndexedFile property for use in rest of code
const { TabixIndexedFile } = window.gmodTABIX

const primateAndPromoterAiTableUrls = {
    '37': 'https://storage.googleapis.com/spliceai-lookup-reference-data/PrimateAI_and_PromoterAI_scores.hg19.20250627.tsv.gz',
    '38': 'https://storage.googleapis.com/spliceai-lookup-reference-data/PrimateAI_and_PromoterAI_scores.hg38.20250627.tsv.gz',
}
const primateAndPromoterAiTables = {
    '37': new TabixIndexedFile({ url: primateAndPromoterAiTableUrls['37'], tbiUrl: `${primateAndPromoterAiTableUrls['37']}.tbi` }),
    '38': new TabixIndexedFile({ url: primateAndPromoterAiTableUrls['38'], tbiUrl: `${primateAndPromoterAiTableUrls['38']}.tbi` }),
}

const alphaMissenseTableUrls = {
    '37': 'https://storage.googleapis.com/spliceai-lookup-reference-data/AlphaMissense_hg19.tsv.gz',
    '38': 'https://storage.googleapis.com/spliceai-lookup-reference-data/AlphaMissense_hg38.tsv.gz',
}
const alphaMissenseTables  = {
    '37': new TabixIndexedFile({ url: alphaMissenseTableUrls['37'], tbiUrl: `${alphaMissenseTableUrls['37']}.tbi` }),
    '38': new TabixIndexedFile({ url: alphaMissenseTableUrls['38'], tbiUrl: `${alphaMissenseTableUrls['38']}.tbi` }),
}

const getUCSCBrowserUrl = (genomeVersion, chrom, pos) => {
    const genomeVers = genomeVersion.replace('37', '19')
    const chromWithoutPrefix = chrom.toUpperCase().replace('CHR', '')

    return `https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg${genomeVers}&position=chr${chromWithoutPrefix}:${pos}`
}

const showError = (errorMessage) => {
    $("#error-box").html( $("#error-box").html() + `<br />${errorMessage.toString()}<br />`)
    $("#error-box").show()
}

const normalizeVariant = async (variant, genomeVersion) => {
    /* Convert the given variant to a standardized "{chrom}-{pos}-{ref}-{alt}" string using a regular expression, or
        * if that fails, assume the variant is in HGVS notation and try using the Ensembl hgvs API to convert it.
        *
        * Args:
        *  variant (string): user input text
        *  genomeVersion (string): "37" or "38"
        *
        * Return:
        *   dictionary with keys:
        *   variant: the input variant reformatted as a "{chrom}-{pos}-{ref}-{alt}" string
        *   consequence: the most severe consequence of the variant according to Ensembl VEP, 
        *   which can be used to prioritize transcripts in the results table
        *  
        */

    const ensemblApiPrefix = `https://${genomeVersion == '37' ? 'grch37.' : ''}rest.ensembl.org/vep/human/hgvs/`

    let chrom, pos, ref, alt, variantConsequence
    let matchedRegExp = variant.match(VARIANT_RE)
    let ensemblApiUrl, ensemblApiResponse, ensemblApiResponseJson
    if (matchedRegExp) {
        // was able to parse the user input using a simple reg-exp
        chrom = matchedRegExp[2].toUpperCase()
        pos = parseInt(matchedRegExp[3].replace(/,/g, ""))
        ref = matchedRegExp[4].toUpperCase()
        alt = matchedRegExp[5].toUpperCase()

        if (ref.length == 1) {
            //handle SNPs and insertions
            ensemblApiUrl = `${ensemblApiPrefix}${chrom}:g.${pos}${ref}>${alt}`
        } else if(alt.length == 1) {
            //handle deletions
            ensemblApiUrl = `${ensemblApiPrefix}${chrom}:g.${pos+1}_${pos+ref.length-1}del${ref.slice(1)}`
        } else {
            //converting MNP into HGVS is not straight forward, so just return the variant
            return { 'variant': `${chrom}-${pos}-${ref}-${alt}` }
        }
    }
    else {
        //assume the variant is already in HGVS notation
        ensemblApiUrl = `${ensemblApiPrefix}${variant.trim()}`
    }


    // try calling the Ensembl API on the user input
    try {
        try {
            ensemblApiResponse = await makeRequest(ensemblApiUrl + "?content-type=application/json&vcf_string=1")
        } catch (e) {
            console.error(e)
            throw Error(`Ensembl API call failed: Unable to reach server`)
        }

        ensemblApiResponseJson = await ensemblApiResponse
        console.log(`Ensembl API response:`, ensemblApiResponseJson)

        if (!ensemblApiResponse.ok || ensemblApiResponseJson.error) {
            let errorText = `${ensemblApiResponseJson.error}`
            //const unableToParseHGVSErrorMatch = errorText.match("Unable to parse HGVS notation")
            //if (unableToParseHGVSErrorMatch) {
            //    errorText = `Ensembl API is unable to parse the variant ${variant}: ${errorText}`
            //}
            const refAlleleErrorMatch = errorText.match(
                new RegExp("[(]([ACGTRYSWKMBDHVN]+)[)] does not match reference allele given by HGVS notation"))
            if (refAlleleErrorMatch) {
                errorText = `${variant} has an unexpected reference allele. The hg${genomeVersion} reference allele should be ${refAlleleErrorMatch[1]}`
            }

            throw new Error(errorText);
        }

        if (!ensemblApiResponseJson[0] || !ensemblApiResponseJson[0].vcf_string) {
            throw new Error(`Unexpected response: ${ensemblApiResponseJson}`);
        }

        variant = ensemblApiResponseJson[0].vcf_string
        variantConsequence = ensemblApiResponseJson[0].most_severe_consequence
        matchedRegExp = variant.match(VARIANT_RE)
        if (!matchedRegExp) {
            throw new Error(`Unexpected response: ${ensemblApiResponseJson}`)
        }

        chrom = matchedRegExp[2].toUpperCase()
        pos = parseInt(matchedRegExp[3])
        ref = matchedRegExp[4].toUpperCase()
        alt = matchedRegExp[5].toUpperCase()

        const result = {
            'variant': `${chrom}-${pos}-${ref}-${alt}`,
            'consequence': variantConsequence,
        }
        console.log("EnsemblAPI result:", result)

        return result

    } catch (e) {
        console.error(e)
        throw new Error(e)
    }
}

const makeRequest = (url) => {
    const method = "GET"
    return new Promise(async (resolve, reject) => {
        const xhr = new XMLHttpRequest()
        xhr.open(method, url)
        xhr.addEventListener("load", () => {
            console.log("GET", url, xhr.status, xhr.statusText)
            let response = {}
            try {
                response = JSON.parse(xhr.response)
            } catch(e) {
                console.error("Unable to parse response", xhr.response)
                reject(`Unexpected error: ${url}`)
                return
            }

            response.ok = xhr.status >= 200 && xhr.status < 300
            resolve(response)
        })
        xhr.addEventListener("error", () => {
            if (xhr.status == 0) {
                reject("Unable to reach server")
            } else {
                reject(`${xhr.status}  ${xhr.statusText}`)
            }
        })
        //console.log('Sending request', method, url)
        xhr.send()
    })
}

const considerInsertedBases = (score, position, ref, alt, scoresForInsertedBases) => {
    //check if the variant is an insertion and if so, return the number of inserted bases
    return score >= 0.01 && position == 0 && ref.length == 1 && alt.length > 1 && scoresForInsertedBases && scoresForInsertedBases.length > 0
}

// TODO: Is it used?
// Special function in case for inserted bases...
const generateTableOfScoresForInsertedBases = (modalId, score, position, ref, alt, scoresForInsertedBases) => {
    //warn about insertion variants that are predicted to cause a donor or acceptor gain somewhere within the inserted sequence.
    if (!considerInsertedBases(score, position, ref, alt, scoresForInsertedBases)) {
        return ""
    }

    let tableRows = []
    for (const scoreObj of scoresForInsertedBases) {
        const RA = parseFloat(scoreObj.RA)
        const RD = parseFloat(scoreObj.RD)
        const AA = parseFloat(scoreObj.AA)
        const AD = parseFloat(scoreObj.AD)

        tableRows.push(
            `<tr><td>${scoreObj.chrom}</td>
                    <td style='text-align:right'>${scoreObj.pos}</td>
                    <td>${scoreObj.ref}</td>
                    <td>${scoreObj.alt}</td>
                    <td ${getScoreStyle(RA)}>${RA}</td>
                    <td ${getScoreStyle(RD)}>${RD}</td>
                    <td ${getScoreStyle(AA)}>${AA}</td>
                    <td ${getScoreStyle(AD)}>${AD}</td></tr>`)
    }

    const table = (
        `<div>
            Detailed SpliceAI predictions for bases within and around the inserted sequence: <br>
            <table class='ui celled table'>
            <thead>
                <tr style='text-align:center'>
                    <th style='width:1%'>chrom</th>
                    <th>position</th>
                    <th style='width:1%'>REF</th>
                    <th style='width:1%'>ALT</th>
                    <th style='width:1%'>REF acceptor score</th>
                    <th style='width:1%'>REF donor score</th>
                    <th style='width:1%'>ALT acceptor score</th>
                    <th style='width:1%'>ALT donor score</th>
                </tr>
            </thead>
            <tbody>
                ${tableRows.join("\n")}
            </tbody>
        </table></div>`
    ).replace(/\n/g, "")

    //check if browser is chrome, safari or firefox
    //const isChrome = /Chrome/.test(navigator.userAgent) && /Google Inc/.test(navigator.vendor)
    //const isSafari = /Safari/.test(navigator.userAgent) && /Apple Computer/.test(navigator.vendor)
    //const isFirefox = /Firefox/.test(navigator.userAgent)

    const popupIcon = `<i style="margin-left:10px; color:#000000" class="table icon score-table" data-position="right center" data-html="${table}"></i>`
    const modal = `<div id='modal-with-table${modalId}' class='ui modal'><i class='close icon'></i><div class='scrolling content'>${table}</div></div>`
    const modalIcon = `<i id='table-icon${modalId}' style='margin-left:10px;color:#000080;cursor:pointer' class='window maximize outline icon'></i>`
    return `${modalIcon}${modal} ${popupIcon}`

}

const updatePositionAccountingForInsertedBases = (scoreKey, score, position, ref, alt, scoresForInsertedBases) => {
    if (!considerInsertedBases(score, position, ref, alt, scoresForInsertedBases)) {
        return `${position} bp`
    }

    if (scoreKey == 'DS_AG') {
        scoreKey = 'AA'
    } else if (scoreKey == 'DS_DG') {
        scoreKey = 'AD'
    } else {
        return `${position} bp`
    }
    // find the position within scoresForInsertedBases at which the score is the highest. Look at the AA score if scoreKey == 'DS_AG', and AD score if scoreKey == 'DS_DG'
    let maxScore = -1
    let maxScorePosition = 0
    for (const scoreObj of scoresForInsertedBases) {
        const scoreValue = parseFloat(scoreObj[scoreKey])
        if (scoreValue > maxScore) {
            maxScore = scoreValue
            maxScorePosition = parseInt(scoreObj.pos)
        }
    }

    return `+${maxScorePosition} bp position within the inserted sequence`
}

const getGnomadDataVersion = (genomeVersion) => {
    return genomeVersion === "38" ? "gnomad_r4" : "gnomad_r2_1"
}

const fetchSplicingToolJson = async (normalizedVariant, variantConsequence, tool, variant, genomeVersion, basicOrComprehensive, maxDistance, mask) => {
    const variantTokens = (normalizedVariant || "---").split("-")
    const chrom = variantTokens[0]
    const pos = variantTokens[1]
    const ref = variantTokens[2]
    const alt = variantTokens[3]

    const baseUrl = baseApiUrl[`${tool.toLowerCase()}-${genomeVersion}`]
    const urlArgs = `hg=${genomeVersion}&bc=${basicOrComprehensive}&distance=${maxDistance}&mask=${mask}&variant=${chrom}-${pos}-${ref}-${alt}&raw=${variant}&variant_consequence=${variantConsequence}`

    let apiResponse
    try {
        apiResponse = await makeRequest(`${baseUrl}/${tool.toLowerCase()}/?${urlArgs}`)
    } catch(e) {
        throw Error(`${tool} API call failed: ${e}`)
    }

    const apiResponseJson = await apiResponse
    console.log(`${tool} API response:`, apiResponseJson)

    if (!apiResponse.ok || apiResponseJson.error) {
        throw Error(`${tool} API call error ${apiResponseJson.error ? `: ${apiResponseJson.error}` : ""}`)
    }

    return apiResponseJson
}

const renderSplicingResultsFromApiJson = (apiResponseJson, normalizedVariant, variantConsequence, tool, variant, genomeVersion, basicOrComprehensive, maxDistance, mask, showRefAltScoreColumns) => {
    /* Render SpliceAI or Pangolin table from an API response JSON (used for single lookups and batch navigation). */

    variantConsequence = variantConsequence || ""
    const variantConsequenceDiv = variantConsequence ? `
        <div display: inline-block>
            <a href="https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html" class="small-link" target="_blank">
            ${variantConsequence.replace(/_/g, ' ')}
            </a>
            </div>
            <br class="only-large-screen"/>
            ` : ""

    //sort transcripts
    apiResponseJson.scores = _(apiResponseJson.scores).sortBy((s) => (  //compute sort key
            100*(s['t_priority'].startsWith("M") ? 0 : 1) +
            10* (s['t_type'] == "protein_coding" ? 0 : 1) +
            -1*  TRANSCRIPT_PRIORITY[s['t_priority']]
    )).values()


    const TRANSCRIPT_PRIORITY_VIEW_LOOKUP = {
        "MS": `<a href="https://www.ncbi.nlm.nih.gov/refseq/MANE" target="_blank">MANE Select transcript</a>`,
        "MP": `<a href="https://www.ncbi.nlm.nih.gov/refseq/MANE" target="_blank">MANE Plus Clinical transcript</a>`,
        "C": "Canonical transcript",
    }

    const noDifferenceDueToNormalization = variant.toLowerCase().trim().replace(/^chr/, "").replace(/[>: _-]+/g, "-") == normalizedVariant.toLowerCase().trim().replace(/^chr/, "").replace(/[>: _-]+/g, "-")
    const normalizedVariantDiv = noDifferenceDueToNormalization ? "" :
        `<br class="only-large-screen"/>
            <div style="margin-left:5px;margin-right:10px; display: inline-block; color:#333333">
            <i>⇒ ${chrom}:${pos} ${ref}&gt;${alt}</i>
            </div>
            <br class="only-large-screen"/>`


    const transcriptCategories = {}
    const tableRows = []
    const modalDialogIds = []
    let transcriptIndex = 0

    $(`#${tool.toLowerCase()}-header`).nextAll().remove()

    for (const scores of apiResponseJson.scores) {
        const isPangolin = tool.toLowerCase() == "pangolin"
        const subRowCount = isPangolin ? 2 : 4
        const strand = scores['t_strand'] == "-" ? "minus" : "plus"
        const gnomadVersion = getGnomadDataVersion(genomeVersion)
        const isMainTranscript = scores['t_priority'] != "N" && (scores['t_priority'] != "C" || transcriptCategories["MS"] == undefined)
        const refSeqLink = scores['t_refseq_ids']? `/ <a href="https://www.ncbi.nlm.nih.gov/search/all/?term=${scores['t_refseq_ids'][0]}" target="_blank">${scores['t_refseq_ids'][0]}</a>` : ""

        transcriptCategories[scores['t_priority']] = true

        const resultRowClasses = [`${tool.toLowerCase()}-result-row`]
        if (isMainTranscript) {
            resultRowClasses.push("main-transcript")
        } else {
            resultRowClasses.push("non-main-transcript")
        }
        if(scores['t_type'] == "protein_coding") {
            resultRowClasses.push("coding-transcript")
        } else {
            resultRowClasses.push("non-coding-transcript")
        }

        let row = `<tr class="${resultRowClasses.join(' ')}">
            <td rowspan="${subRowCount}" style="vertical-align: top">
                <div style="margin-right:10px; display: inline-block">${variant}</div><br class="only-large-screen"/>
                ${normalizedVariantDiv}
                <br class="only-large-screen"/>
                ${variantConsequenceDiv}
                <a href="${getUCSCBrowserUrl(genomeVersion, chrom, pos)}" class="small-link" target="_blank">UCSC</a>,
                <a href="https://gnomad.broadinstitute.org/variant/${chrom}-${pos}-${ref}-${alt}?dataset=${gnomadVersion}" class="small-link" target="_blank">gnomAD</a>
                <br class="only-large-screen"/>

            </td>
            <td rowspan="${subRowCount}" style="vertical-align: top">
                <div style="margin-right:10px; display: inline-block">
                    ${scores['g_name']}
                        <div class="small-link"> (
                            <a href="https://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=${scores['g_id'].split('.')[0]}" target="_blank">${scores['g_id']}</a>
                                / <a href="https://useast.ensembl.org/Homo_sapiens/Transcript/Summary?t=${scores['t_id'].split('.')[0]}" target="_blank">${scores['t_id']}</a>
                                ${refSeqLink})
                        </div><br />
                        <br class="only-large-screen" />
                    <div class="small-link"><a href="https://www.gencodegenes.org/pages/biotypes.html" target="_blank">${scores['t_type'].replace(/_/g, " ")}</a></div>
                    <div class="small-link">${scores['t_priority'] != "N" ? TRANSCRIPT_PRIORITY_VIEW_LOOKUP[scores['t_priority']] : ""}</div>
                    <div class="small-link"> (${strand} strand)</div>
                </div><br class="only-large-screen" />
                <br class="only-large-screen" />
                <a href="https://www.omim.org/search?search=${scores['g_name']}" class="small-link" target="_blank">OMIM</a>,
                <a href="https://gtexportal.org/home/gene/${scores['g_name']}" class="small-link" target="_blank">GTEx</a>,
                <a href="https://gnomad.broadinstitute.org/gene/${scores['g_name']}?dataset=${gnomadVersion}" class="small-link" target="_blank">gnomAD</a>,
                <a href="https://search.clinicalgenome.org/kb/genes?page=1&size=25&search=${scores['g_name']}" class="small-link" target="_blank">ClinGen</a>,
                <a href="https://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=${scores['g_name']}" class="small-link" target="_blank">Ensembl</a>,
                <a href="https://www.deciphergenomics.org/gene/${scores['g_name']}" class="small-link" target="_blank">Decipher</a>,
                <a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=${scores['g_name']}" class="small-link" target="_blank">GeneCards</a>
            </td>`

        for (const [i, label, scoreKey, positionKey, REF_scoreKey, ALT_scoreKey] of (
            isPangolin ? [
                [0, "Splice Loss", "DS_SL", "DP_SL", "SL_REF", "SL_ALT"],
                [1, "Splice Gain", "DS_SG", "DP_SG", "SG_REF", "SG_ALT"],
            ] : [
                [0, "Acceptor Loss", "DS_AL", "DP_AL", "DS_AL_REF", "DS_AL_ALT"],
                [1, "Donor Loss",    "DS_DL", "DP_DL", "DS_DL_REF", "DS_DL_ALT"],
                [2, "Acceptor Gain", "DS_AG", "DP_AG", "DS_AG_REF", "DS_AG_ALT"],
                [3, "Donor Gain",    "DS_DG", "DP_DG", "DS_DG_REF", "DS_DG_ALT"],
            ])) {
            if (i > 0) {
                row += `</tr><tr class="${resultRowClasses.join(' ')}">`
            }

            const modalDialogId = transcriptIndex*10 + i
            const detailedScoresTable = generateTableOfScoresForInsertedBases(
                modalDialogId, scores[scoreKey], scores[positionKey], ref, alt, scores["SCORES_FOR_INSERTED_BASES"])

            if (detailedScoresTable) {
                modalDialogIds.push(modalDialogId)
            }

            row += `<td>${label}</td>
                    <td ${getScoreStyle(scores[scoreKey])}>
                        ${formatScore(scores[scoreKey])}
                        ${detailedScoresTable}
                        </td>
                    <td>${
                        (
                            parseFloat(scores[scoreKey]) == 0 &&
                            parseFloat(scores[REF_scoreKey]) == 0 &&
                            parseFloat(scores[ALT_scoreKey]) == 0
                        )? "" : updatePositionAccountingForInsertedBases(scoreKey, scores[scoreKey], scores[positionKey], ref, alt, scores["SCORES_FOR_INSERTED_BASES"])}
                    </td>`

            if (showRefAltScoreColumns == "1") {
                row += `<td class="ref-score-column">${formatScore(scores[REF_scoreKey])}</td>
                        <td class="alt-score-column">${formatScore(scores[ALT_scoreKey])}</td>`
            }
        }
        row += "</tr>"

        tableRows.push(row)
        transcriptIndex++
    }

    $(`#${tool.toLowerCase()}-header`).after(tableRows.join(""))

    if (transcriptCategories["MS"]) {
        $(".main-transcript-label").html("MANE Select")
    } else if (transcriptCategories["MP"]) {
        $(".main-transcript-label").html("MANE Plus Clinical")
    } else if (transcriptCategories["C"]) {
        $(".main-transcript-label").html("Canonical")
    } else {
        $("#transcript-button-table").hide()
    }

    // initialize any modal dialogs
    for (const modalDialogIndex of modalDialogIds) {
        $(`#table-icon${modalDialogIndex}`).click(() => {
            $(`#modal-with-table${modalDialogIndex}`).modal('show')
        })
    }
    $("#transcript-button-table").show()
    updateTranscriptButtons(Object.keys(transcriptCategories).length > 1 ? "main" : "all")
}

const generateSplicingResultsTable = async (normalizedVariant, variantConsequence, tool, variant, genomeVersion, basicOrComprehensive, maxDistance, mask, showRefAltScoreColumns) => {
    const apiResponseJson = await fetchSplicingToolJson(normalizedVariant, variantConsequence, tool, variant, genomeVersion, basicOrComprehensive, maxDistance, mask)
    renderSplicingResultsFromApiJson(apiResponseJson, normalizedVariant, variantConsequence, tool, variant, genomeVersion, basicOrComprehensive, maxDistance, mask, showRefAltScoreColumns)
    return apiResponseJson
}


const generateOtherPredictorsTable = async (normalizedVariant, variantConsequence, variant, genomeVersion) => {
    /* Generate the results table to show either the SpliceAI or Pangolin scores */
    console.log(`Generating other scores table for ${normalizedVariant} from ${primateAndPromoterAiTableUrls[genomeVersion]}`)

    const variantTokens = (normalizedVariant || "---").split("-")
    const chrom = `chr${variantTokens[0].replace("chr", "")}`
    const pos = parseInt(variantTokens[1])
    const ref = variantTokens[2]
    const alt = variantTokens[3]

    //example variant: 1-55039916-G-A  (hg38)

    const otherPredictorScores = {}
    const queryLookupTable = async () => {
        await Promise.all([
            primateAndPromoterAiTables[genomeVersion].getLines(chrom, pos-1, pos, line => {
                //console.log("Got line:", line)
                const fields = line.split('\t')
                const lineChrom = fields[0]
                const linePos = parseInt(fields[1])
                const lineRef = fields[2]
                const lineAlt = fields[3]
                //console.log("`Got scores`:", fields)

                if (linePos != pos || lineRef != ref || lineAlt != alt) {
                    return
                }
                const percentile = parseFloat(fields[4])
                const genePercentilethreshold = parseFloat(fields[5])
                //console.log("PrimateAI-3D", percentile, genePercentilethreshold)
                if (!isNaN(percentile)) {
                    otherPredictorScores['primateai3d']  = {
                        'percentile': percentile,
                        'genePercentileThreshold': genePercentilethreshold,
                    }
                }
                const promoterAiScore = parseFloat(fields[6])
                //console.log("PromoterAI", promoterAiScore)
                if (!isNaN(promoterAiScore)) {
                    otherPredictorScores['promoterai'] = {
                        'score': promoterAiScore,
                    }
                }
            },),

            alphaMissenseTables[genomeVersion].getLines(chrom, pos-1, pos, line => {
                //console.log("Got line:", line)
                const fields = line.split('\t')
                const lineChrom = fields[0]
                const linePos = parseInt(fields[1])
                const lineRef = fields[2]
                const lineAlt = fields[3]
                //console.log("`Got scores`:", fields)

                if (linePos != pos || lineRef != ref || lineAlt != alt) {
                    return
                }
                //const transcriptId = fields[6]
                //const proteinVariant = fields[7]
                //const consequence = fields[9]  //benign, ambiguous, pathogenic
                const alphaMissenseScore = parseFloat(fields[8])

                //console.log("AlphaMissense", alphaMissenseScore)
                if (!isNaN(alphaMissenseScore)) {
                    otherPredictorScores['alphamissense'] = {
                        'score': alphaMissenseScore,
                        //'transcriptId': transcriptId,
                        //'proteinVariant': proteinVariant,
                        //'consequence': consequence,
                    }
                }
            },)
        ])
    }

    const queryGnomAD = async () => {
        console.log("Querying gnomAD")
        query = `{
            variant(variantId: "${chrom.replace("chr", "")}-${pos}-${ref}-${alt}", dataset: ${getGnomadDataVersion(genomeVersion)}) {
            in_silico_predictors {
                id,
                value
            },
            }
        }`

        let response
        try {
            response = await fetch("https://gnomad.broadinstitute.org/api", {
                "method": "POST",
                "headers": {
                    "Content-Type": "application/json",
                    "Accept": "application/json",
                },
                "body": JSON.stringify({ query }),
            }) //, { mode: 'no-cors' })
        } catch (e) {
            console.error("gnomAD query failed:", e)
            return
        }


        if (response.ok) {
            const responseJson = await response.json()
            console.log("gnomAD response:", responseJson)
            if (responseJson && responseJson.data && responseJson.data.variant) {
                /*
                if (responseJson.data.variant.joint && responseJson.data.variant.joint.ac && responseJson.data.variant.joint.an) {
                    otherPredictorScores.push({
                        'name': 'gnomAD',
                        'score': `AC: ${responseJson.data.variant.joint.ac}, &nbsp; AN: ${responseJson.data.variant.joint.an}`,
                    })
                }
                */

                if (responseJson.data.variant.in_silico_predictors) {
                    responseJson.data.variant.in_silico_predictors.forEach((predictor) => {
                        if (predictor.id.startsWith("spliceai") || predictor.id.startsWith("pangolin")) {
                            return
                        }
                        const predictorName = predictor.id
                        const score = parseFloat(parseFloat(predictor.value).toFixed(3))
                        if (otherPredictorScores[predictorName] && Math.abs(otherPredictorScores[predictorName].score - score) > 0.001) {
                            console.warn(`Mismatch between gnomAD and myvariant.info scores for ${predictorName}: ${score} vs ${otherPredictorScores[predictorName].score}`)
                        }
                        otherPredictorScores[predictorName] = {
                            'score': score,
                        }
                    })
                }
            }
        }
    }

    const queryMyVariantInfo = async () => {
        console.log("Querying myvariant.info")
        const variantId = `${chrom}:g.${pos}${ref}>${alt}`
        //add a data section to the body of the request
        let response
        try {
            response = await fetch(`https://myvariant.info/v1/variant/${encodeURIComponent(variantId)}?fields=dbnsfp&size=10&assembly=${genomeVersion == '37' ? 'hg19' : 'hg38'}`, {
                "method": "GET",
                "headers": {
                    "accept": "*/*",
                }
            })
        } catch (e) {
            console.error("myvariant.info query failed:", e)
            return
        }

        if (response.ok) {
            const responseJson = await response.json()
            console.log("myvariant.info response:", responseJson)
            if (responseJson.dbnsfp) {
                for (const otherPredictorName of [
                    "revel.score", "cadd.phred", "alphamissense.score", "sift.score",
                ]) {
                    const keyTokens = otherPredictorName.split(".")
                    const key1 = keyTokens[0]
                    const key2 = keyTokens[1]
                    if (responseJson.dbnsfp[key1] && responseJson.dbnsfp[key1][key2] != null) {
                        let score = responseJson.dbnsfp[key1][key2]
                        if (Array.isArray(score)) {
                            score = score[0]  // if its an array, get the first elegment
                        }
                        const otherPredictorName = key1.replace("revel", "revel_max").replace("sift", "sift_max")
                        if (otherPredictorScores[otherPredictorName]) {
                            if (Math.abs(otherPredictorScores[otherPredictorName].score - score) > 0.001) {
                                console.error(`Mismatch between myvariant.info and gnomAD scores for ${key1}: ${score} vs ${otherPredictorScores[otherPredictorName].score}`)
                            }
                        } else {
                            otherPredictorScores[otherPredictorName] = {
                                'score': parseFloat(parseFloat(score).toFixed(3)),
                            }
                        }
                    }
                }
            }
        }
    }

    console.log("Querying lookup tables, gnomAD, and myvariant.info")
    await Promise.all([queryLookupTable(), queryGnomAD(), queryMyVariantInfo()]) //.map(p => p.catch(e => ({'error': e.message})))

    console.log("Other predictors scores:", otherPredictorScores)
    const noDifferenceDueToNormalization = variant.toLowerCase().trim().replace(/^chr/, "").replace(/[>: _-]+/g, "-") == normalizedVariant.toLowerCase().trim().replace(/^chr/, "").replace(/[>: _-]+/g, "-")
    const normalizedVariantDiv = noDifferenceDueToNormalization ? "" :
        `<br class="only-large-screen"/>
            <div style="margin-left:5px;margin-right:10px; display: inline-block; color:#333333">
            <i>⇒ ${chrom}:${pos} ${ref}&gt;${alt}</i>
            </div>
            <br class="only-large-screen"/>`


    const otherPredictorNames = [
        "alphamissense", "cadd", "phylop", "polyphen_max", "primateai3d", "promoterai", "revel_max", "sift_max",
    ]

    const otherPredictorMap = {}
    for (const [otherPredictorName, otherPredictor] of Object.entries(otherPredictorScores)) {
        otherPredictorMap[otherPredictorName] = otherPredictor
    }

    $(`#other-predictors-header`).nextAll().remove()
    let firstColumn = `
        <td class="eight wide column" style="vertical-align: top" rowSpan="${Object.keys(otherPredictorScores).length}">
            <div style="margin-right:10px; display: inline-block">${variant}</div><br class="only-large-screen"/>
            ${normalizedVariantDiv}
        </td>`


    const tableRows = []
    const missingScores = []
    for (const otherPredictorName of otherPredictorNames) {
        const predictorLabel = nameMapForPredictorScores[otherPredictorName] || otherPredictorName
        const otherPredictor = otherPredictorMap[otherPredictorName]
        const bgColor = otherPredictor && colorMapForPredictorScores[otherPredictorName] ? colorMapForPredictorScores[otherPredictorName](otherPredictor) : "#ffffff"
        const helpIcon = otherPredictor && helpTextForPredictorScores[otherPredictorName] ? `<i class='question circle outline icon' data-html="${helpTextForPredictorScores[otherPredictorName](otherPredictor)}"></i>` : ""

        let row = `<tr>${firstColumn}<td class="four wide column" style="vertical-align: top; white-space: nowrap">${predictorLabel}</td>`
        if (otherPredictorName == "primateai3d") {
            if (otherPredictor) {
                row += `<td class="three wide column" style="vertical-align: top; white-space: nowrap; background-color: ${bgColor}">${otherPredictor.percentile.toFixed(2)} &nbsp; (gene-specific threshold: ${otherPredictor.genePercentileThreshold.toFixed(2)}) &nbsp; ${helpIcon}</td><td class="one wide column" style="background-color: ${bgColor}"></td></tr>`
                tableRows.push(row)
                firstColumn = ""
            } else {
                missingScores.push(otherPredictorName)
            }
        } else if (otherPredictorName == "promoterai") {
            if (otherPredictor) {
                row += `<td class="three wide column" style="vertical-align: top; white-space: nowrap; background-color: ${bgColor}">${otherPredictor.score.toFixed(2)} &nbsp; ${helpIcon}</td>
                <td class="one wide column" style="vertical-align: top; white-space: nowrap; background-color: ${bgColor}"></td></tr>`
                tableRows.push(row)
                firstColumn = ""
            } else {
                missingScores.push(otherPredictorName)
            }
        } else if (otherPredictorName in predictorScoreToPoints) {
            if (otherPredictor) {
                row += `<td class="three wide column" style="vertical-align: top; white-space: nowrap; background-color: ${bgColor}">${otherPredictor.score}</td>`
                if (variantConsequence && variantConsequence.toLowerCase().includes("missense")) {
                    row += `<td class="one wide column" style="vertical-align: top; text-align: right; white-space: nowrap; background-color: ${bgColor}">
                            ${predictorScoreToPoints[otherPredictorName](otherPredictor)}&nbsp;${helpIcon}
                        </td></tr>`
                } else {
                    row += `<td class="one wide column" style="background-color: ${bgColor}"></td>`
                }
                tableRows.push(row)
                firstColumn = ""
            } else {
                missingScores.push(otherPredictorName)
            }
        } else {
            console.error("Unknown predictor", otherPredictorName)
        }
    }

    if (missingScores.length > 0) {
        console.log("Missing scores for", missingScores)
        const missingScoresRow = `<tr><td class="four wide column"></td><td class="four wide column" style="vertical-align: top; white-space: nowrap" colspan="3"><div style="padding-top: 10px; display:inline-block">${missingScores.map((name, i) => `${i === 0 || i < missingScores.length - 1 ? nameMapForPredictorScores[name] : 'and ' + nameMapForPredictorScores[name]}`).join(", &nbsp;")} scores are not available for this variant &nbsp; <i class='question circle outline icon' data-position='right center' data-content='AlphaMissense, PrimateAI-3D, and PromoterAI scores are retrieved from public lookup tables of precomputed scores, while CADD, PhyloP, PolyPhen, REVEL, and SIFT scores are retrieved from the gnomAD and myvariant.info APIs'/></div></td></tr>`
        tableRows.push(missingScoresRow)
    }
    $("#other-predictors-header").after(tableRows.join(""))
}

const updateTranscriptButtons = (category) => {
    $("#main-transcript-button, #all-transcript-button").removeClass("primary")
    $(`#${category}-transcript-button`).addClass("primary")

    if (category == "main") {
        $(".spliceai-result-row").hide()
        $(".spliceai-result-row.main-transcript").show()
    } else if (category == "all") {
        $(".spliceai-result-row").show()
    }
}

const updateVisualizationCheckboxes = (readFromLocalStorage, genomeVersion) => {
    /** 
     * Retrieve and save the current state of checkboxes in the Visualization section, or set the state of these 
     * checkboxes based on cookies / local storage.
     */
    let tracksToShow
    if (readFromLocalStorage) {
        // default tracks
        tracksToShow = {
            "igv-variant": true,
            "igv-spliceai-ref-alt": true,
            "igv-spliceai-delta-scores": true,
            "igv-pangolin-delta-scores": true,
        }

        // if there are tracksToShow in local storage, overwrite defaults with those settings
        const fromStorage = localStorage.getItem("tracksToShow")
        if (fromStorage) {
            const fromStorageJson = JSON.parse(fromStorage)
            for (const key of Object.keys(fromStorageJson)) {
                tracksToShow[key] = fromStorageJson[key]
                try {
                    $(`input[name='${key}']`).prop( "checked", tracksToShow[key] )
                } catch (e) {
                    console.log(e)
                }
            }
        }
    } else {

        // retrieve state from checkboxes
        tracksToShow = {}
        for (const name of ["igv-variant", "igv-spliceai-ref-alt", "igv-spliceai-delta-scores", "igv-pangolin-delta-scores", "igv-gencode-genes"]) {
            tracksToShow[name] = $(`input[name='${name}']`).prop("checked")
        }
        for (const name of ["igv-100-mer-mappability", "igv-segdups", "igv-mane-genes"]) {
            // these tracks are not available for hg37
            $(`input[name='${name}']`).prop("checked", $(`input[name='${name}']`).prop("checked") && genomeVersion != "37")
            $(`input[name='${name}']`).prop("disabled", genomeVersion == "37")

            tracksToShow[name] = $(`input[name='${name}']`).prop("checked")
        }
        for (const minScore of [0.5, 0.2]) {   //0.2,  hide the >= 0.2 tracks for now
            for (const splicePredictionType of ["loss", "gain"]) {
                const name = `igv-spliceai-precomputed-${splicePredictionType}-${minScore}`

                // these tracks are not available for hg19
                $(`input[name='${name}']`).prop("checked", $(`input[name='${name}']`).prop("checked") && genomeVersion != "37")
                $(`input[name='${name}']`).prop("disabled", genomeVersion == "37")

                tracksToShow[name] = $(`input[name='${name}']`).prop("checked")
            }
        }

        tracksToShow[`igv-spliceai-precomputed-score-genes`] = $(`input[name='igv-spliceai-precomputed-score-genes']`).prop("checked")

        for (const tissue of ["blood", "fibroblasts", "muscle", "lymphocytes", "brain-cortex"]) {
            // these tracks are not available for hg19
            $(`input[name='igv-gtex-${tissue}']`).prop("checked", $(`input[name='igv-gtex-${tissue}']`).prop("checked") && genomeVersion != "37")
            $(`input[name='igv-gtex-${tissue}']`).prop("disabled", genomeVersion == "37")
            $(`input[name='igv-gtex-${tissue}-all']`).prop("checked", $(`input[name='igv-gtex-${tissue}-all']`).prop("checked") && genomeVersion != "37")
            $(`input[name='igv-gtex-${tissue}-all']`).prop("disabled", genomeVersion == "37")

            tracksToShow[`igv-gtex-${tissue}`] = $(`input[name='igv-gtex-${tissue}']`).prop("checked")
            tracksToShow[`igv-gtex-${tissue}-all`] = $(`input[name='igv-gtex-${tissue}-all']`).prop("checked")
        }
        //console.log("Set local storage tracksToShow", tracksToShow)
        localStorage.setItem("tracksToShow", JSON.stringify(tracksToShow))
    }
    
    return tracksToShow
}


const generateIgvConfig = (spliceaiResponseJson, pangolinResponseJson, genomeVersion) => {
    /* Generates an igv.js config json objects based on the predicted scores from SpliceAI and/or Pangolin
    *
    * Args:
    *   allNonZeroScoresFromSpliceAI (array): an array of acceptor/donor loss/gain scores from SpliceAI
    *   allNonZeroScoresFromPangolin (array): an array of acceptor/donor loss/gain scores from Pangolin
    *   genomeVersion (string): "37" or "38"
    *
    * Return:
    *   object: an igv.js config json object
    */

    const apiResponseJson = spliceaiResponseJson || pangolinResponseJson
    let chrom = apiResponseJson.chrom
    chrom = `chr${chrom.replace('chr', '')}`
    const variantPos = apiResponseJson.pos
    const variantRef = apiResponseJson.ref
    const variantAlt = apiResponseJson.alt

    const tracksToShow = updateVisualizationCheckboxes(false, genomeVersion)

    // log event
    makeRequest(`${baseApiUrl['pangolin-37']}/log/show_igv?details=${encodeURIComponent(JSON.stringify(tracksToShow))}&hg=${genomeVersion}&distance=${apiResponseJson.distance}&mask=${apiResponseJson.mask}&variant=${apiResponseJson.variant}`)

    let minPos = 1e9
    let maxPos = 0
    for (const apiResponse of [spliceaiResponseJson, pangolinResponseJson]) {
        if (!apiResponse) {
            continue
        }
        if (!apiResponse.allNonZeroScores) {
            continue
        }
        for (const scores of apiResponse.allNonZeroScores) {
            scores["chr"] = chrom
            scores["start"] = scores["pos"]
            scores["end"] = scores["pos"]
            minPos = Math.min(minPos, scores["pos"])
            maxPos = Math.max(maxPos, scores["pos"])
        }
    }

    // specify IGV tracks and the data to display in them
    const tracks = []

    tracks.push({
        name: "Refseq",
        format: "refgene",
        url: genomeVersion == "38" ?
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz" :
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ncbiRefSeq.txt.gz",
        indexed: false,
        infoURL: "https://www.ncbi.nlm.nih.gov/gene/?term=$$",
        height: 100,
    })

    if (tracksToShow["igv-variant"]) {
        tracks.push({
            name: "Variant",
            type: "vcf",
            description: `This track shows the location of the <b>${chrom}:${variantPos} ${variantRef}>${variantAlt}</b> variant`,
            height: 30,
            features: [
                {
                    chr: chrom,   //see createVCFVariant function in the igv.js repo for the allowed fields
                    pos: variantPos,
                    start: variantPos - 1,
                    end: variantPos + variantRef.length - 1,
                    referenceAllele: variantRef,
                    alternateBases: variantAlt,
                    names: ".", // id in VCF
                    info: {
                        variant: `${chrom}-${variantPos}-${variantRef}-${variantAlt}`,
                    }
                },
            ],
        })
    }

    if (spliceaiResponseJson) {
        if (tracksToShow["igv-spliceai-ref-alt"]) {
            tracks.push({
                name: `SpliceAI REF/ALT`,
                description: `This track shows SpliceAI scores for the
                                <span style="color:#0000B4"><b>reference sequence</b></span> (without the variant) and the
                                <span style="color:#05d0d2"><b>alternate sequence</b></span> (with the variant). <br />
                                An <b>A</b> (acceptor) or <b>D</b> (donor) symbol is shown where SpliceAI predicts there to be a splice donor or acceptor with score ≥ 0.01. <br/>
                                These symbols are shown up-side-down when the score for the alternate sequence is lower than the score for the reference sequence.
                                Numberical labels represent scores for the reference sequence. These predictions are based on the pre-mRNA sequence for transcript <br />
                                <b>${spliceaiResponseJson.allNonZeroScoresTranscriptId}</b> (${spliceaiResponseJson.allNonZeroScoresStrand == '-' ? 'minus' : 'plus'} strand)
                                within a +/- ${spliceaiResponseJson.distance}bp window around <b>${chrom}:${variantPos} ${variantRef}>${variantAlt}</b>.`,
                height: 100,
                rawOrDelta: "raw",
                tool: "SpliceAI",
                type: "spliceprediction",
                features: spliceaiResponseJson.allNonZeroScores || [],
                strand: spliceaiResponseJson.allNonZeroScoresStrand || "+",
            })
        }
        if (tracksToShow["igv-spliceai-delta-scores"]) {
            tracks.push({
                name: `SpliceAI Δ`,
                description: `This track shows <b>A</b> (acceptor) and <b>D</b> (donor) symbols at positions where SpliceAI predicts a delta score ≥ 0.01. <br/>
                                Vertical lines are <b style="color:#FF0000">red</b> for delta scores ≥ 0.8,
                                <b style="color:#FFCB1F">yellow</b> for delta scores ≥ 0.5, and
                                <b style="color:#1fb839">green</b> for delta scores ≥ 0.2<br />
                                These predictions are based on the pre-mRNA sequence for transcript
                                <b>${spliceaiResponseJson.allNonZeroScoresTranscriptId}</b> (${spliceaiResponseJson.allNonZeroScoresStrand == '-' ? 'minus' : 'plus'} strand)<br />
                                within a +/- ${spliceaiResponseJson.distance}bp window around <b>${chrom}:${variantPos} ${variantRef}>${variantAlt}</b>.`,
                height: 200,
                rawOrDelta: "delta",
                tool: "SpliceAI",
                type: "spliceprediction",
                features: spliceaiResponseJson.allNonZeroScores || [],
                strand: spliceaiResponseJson.allNonZeroScoresStrand || "+",
            })
        }
    }

    if (pangolinResponseJson) {
        if (tracksToShow["igv-pangolin-delta-scores"]) {

            tracks.push({
                name: `Pangolin Δ`,
                description: `This track shows a <b>P</b> at positions where Pangolin predicts a delta score ≥ 0.01. <br/>
                                Vertical lines are <b style="color:#FF0000">red</b> for delta scores ≥ 0.8,
                                <b style="color:#FFCB1F">yellow</b> for delta scores ≥ 0.5, and
                                <b style="color:#1fb839">green</b> for delta scores ≥ 0.2<br />
                                These predictions are based on the pre-mRNA sequence for transcript
                                <b>${pangolinResponseJson.allNonZeroScoresTranscriptId}</b> (${pangolinResponseJson.allNonZeroScoresStrand == '-' ? 'minus' : 'plus'} strand)<br />
                                within a +/- ${pangolinResponseJson.distance}bp window around <b>${chrom}:${variantPos} ${variantRef}>${variantAlt}</b>.`,
                height: 200,
                rawOrDelta: "delta",
                tool: "Pangolin",
                type: "spliceprediction",
                features: pangolinResponseJson.allNonZeroScores || [],
                strand: pangolinResponseJson.allNonZeroScoresStrand || "+",
            })
        }
    }

    if (tracksToShow["igv-gencode-genes"]) {
        const gencodeTrackPath = `gs://tgg-viewer/ref/GRCh${genomeVersion}/gencode_${GENCODE_VERSION}/gencode.${GENCODE_VERSION}.GRCh${genomeVersion}.sorted.txt.gz`

        tracks.push({
            name: `Gencode ${GENCODE_VERSION}`,
            format: 'refgene',
            url: gencodeTrackPath,
            indexURL: `${gencodeTrackPath}.tbi`,
            indexed: true,
            searchable: true,
            height: 350,
            visibilityWindow: -1,
            order: 1000001,
            displayMode: 'EXPANDED',
            color: 'rgb(76,171,225)',
        })
    }

    if (tracksToShow["igv-mane-genes"]) {
        tracks.push({
            name: "MANE v1.4",
            format: "gtf",
            url: "gs://tgg-viewer/ref/GRCh38/MANE_v1_4/MANE.GRCh38.v1.4.ensembl_genomic.sorted.gtf.gz",
            indexURL: "gs://tgg-viewer/ref/GRCh38/MANE_v1_4/MANE.GRCh38.v1.4.ensembl_genomic.sorted.gtf.gz.tbi",
            height: 100,
        })
    }

    if (tracksToShow["igv-spliceai-precomputed-score-genes"]) {
        const precomputedScoresGeneTrackPath = `gs://tgg-viewer/ref/GRCh${genomeVersion}/gencode_v24/gencode_v24_annotations.grch${genomeVersion}.bed.gz`

        tracks.push({
            name: `Genes used for precomputed scores`,
            //format: 'bed',
            url: `${precomputedScoresGeneTrackPath}?`,
            indexURL: `${precomputedScoresGeneTrackPath}.tbi?`,
            indexed: true,
            searchable: true,
            height: 350,
            visibilityWindow: -1,
            displayMode: 'EXPANDED',
            color: 'rgb(76,171,225)',
        })
    }


    if (genomeVersion == "38") {
        for (const minScore of [0.5, 0.2]) {
            for (const splicePredictionType of ["loss", "gain"]) {
                if (tracksToShow[`igv-spliceai-precomputed-${splicePredictionType}-${minScore}`]) {
                    tracks.push(
                        {
                            name: `SpliceAI: A or D ${splicePredictionType} ≥ ${minScore}`,
                            description: `This track visualizes Illumina's precomputed SpliceAI score tables for SNVs and small INDELs.<br/>
                                    Each genomic location of an SNV or INDEL variant with a SpliceAI &#916; score ≥ ${minScore} is shown as the origin of an arrow.<br />
                                    The arrow points to the location where that variant would likely cause acceptor or donor ${splicePredictionType}.
                                    Clicking on the arrow displays additional details.`,
                            type: "spliceJunctions",
                            height: 100,
                            url: `gs://tgg-viewer/ref/GRCh38/spliceai/spliceai_scores.raw.snps_and_indels.hg38.filtered.sorted.score_${minScore}.splice_${splicePredictionType}.bed.gz`,
                            indexURL: `gs://tgg-viewer/ref/GRCh38/spliceai/spliceai_scores.raw.snps_and_indels.hg38.filtered.sorted.score_${minScore}.splice_${splicePredictionType}.bed.gz.tbi`,
                        }
                    )
                }
            }
        }


        for (const [filenamePrefix, tissue, sampleCount] of [
            ["GTEX_blood.755_samples", "blood", 755],
            ["GTEX_fibs.504_samples", "fibroblasts", 504],
            ["GTEX_muscle.803_samples", "muscle", 803],
            ["GTEX_lymphocytes.174_samples", "lymphocytes", 174],
            ["GTEX_brain_cortex.255_samples", "brain-cortex", 255],
            //["GTEX_frontal_cortex.209_samples", "frontal cortex", 803],
        ]) {
            for (const normalized of [true, false]) {
                if (!tracksToShow[`igv-gtex-${tissue}${normalized ? '' : '-all'}`]) {
                    continue
                }

                const filenameSuffix = normalized ? ".normalized" : ""
                tracks.push({
                    name: `GTEx ${tissue}: ${normalized ? 'per-sample average' : `summed over all ${sampleCount} samples`}`,
                    description: `This track shows the combined splice junctions from all ${sampleCount} ${tissue} samples available in GTEx v8.
                                    Splice junctions are labeled with the total number of RNA-seq reads that supported the junction, summed across the ${sampleCount} samples${normalized ? ' and divided by ' + sampleCount : ''}.
                                    The coverage track shows the total number of RNA-seq reads that overlap each position, summed across all samples${normalized ? ' and divided by ' + sampleCount : ''}.`,
                    type: "merged",
                    height: 100,
                    tracks: [{
                        type: "wig",
                        format: "bigwig",
                        url: `gs://tgg-viewer/ref/GRCh38/gtex_v8/${filenamePrefix}${filenameSuffix}.bigWig`,
                    }, {
                        type: 'spliceJunctions',
                        format: 'bed',
                        minJunctionEndsVisible: 1,
                        colorBy: 'strand',
                        minTotalReads: 3,
                        url: `gs://tgg-viewer/ref/GRCh38/gtex_v8/${filenamePrefix}${filenameSuffix}.junctions.bed.gz`,
                        indexURL: `gs://tgg-viewer/ref/GRCh38/gtex_v8/${filenamePrefix}${filenameSuffix}.junctions.bed.gz.tbi`,
                    }],
                })
            }
        }

        if (tracksToShow["igv-100-mer-mappability"]) {
            tracks.push({
                name: "100-mer mappability",
                description: "100bp k-mer mappability track from UCSC. Lower values indicate regions that are not unique",
                type: "wig",
                format: "bigwig",
                url: "gs://tgg-viewer/ref/GRCh38/mappability/GRCh38_no_alt_analysis_set_GCA_000001405.15-k100_m2.bw",
                height: 100,
            })
        }

        if (tracksToShow["igv-segdups"]) {
            tracks.push({
                name: "SegDups",
                description: "Segmental duplications track from UCSC",
                format: "gtf",
                url: "gs://tgg-viewer/ref/GRCh38/segdups/segdups.gtf.gz",
                indexURL: "gs://tgg-viewer/ref/GRCh38/segdups/segdups.gtf.gz.tbi",
                height: 100,
            })
        }
    }

    let locusMargin = maxPos - minPos < 200 ? 15 : 50

    // igv.js reference genome specs are copied from https://igv.org/genomes/genomes.json
    // These must be specified explicitly rather than just setting the reference to "hg19" or "hg38" so that
    // the RefSeq gene track can be the first track (as described in https://github.com/igvteam/igv.js/issues/1518)
    let reference = genomeVersion == "37" ? {
        "id": "hg19",
        "name": "Human (GRCh37/hg19)",
        "fastaURL": "https://igv-genepattern-org.s3.amazonaws.com/genomes/seq/hg19/hg19.fasta",
        "indexURL": "https://igv-genepattern-org.s3.amazonaws.com/genomes/seq/hg19/hg19.fasta.fai",
        "cytobandURL": "https://igv-genepattern-org.s3.amazonaws.com/genomes/seq/hg19/cytoBand.txt",
        "aliasURL": "https://s3.amazonaws.com/igv.org.genomes/hg19/hg19_alias.tab",
        "chromosomeOrder": "chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chrX, chrY",
        "tracks": tracks,
    } : {
        "id": "hg38",
        "name": "Human (GRCh38/hg38)",
        "fastaURL": "https://igv-genepattern-org.s3.amazonaws.com/genomes/seq/hg38/hg38.fa",
        "indexURL": "https://igv-genepattern-org.s3.amazonaws.com/genomes/seq/hg38/hg38.fa.fai",
        "cytobandURL": "https://s3.amazonaws.com/igv.org.genomes/hg38/annotations/cytoBandIdeo.txt.gz",
        "aliasURL": "https://s3.amazonaws.com/igv.org.genomes/hg38/hg38_alias.tab",
        "chromosomeOrder": "chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chrX, chrY",
        "tracks": tracks,
    }

    return {
        reference: reference,
        showCursorTrackingGuide: true,
        locus: `${chrom}:${minPos - locusMargin}-${maxPos + locusMargin}`,
    }
}

const readFileAsText = (file) => new Promise((resolve, reject) => {
    const r = new FileReader()
    r.onload = () => resolve(r.result)
    r.onerror = () => reject(new Error("Could not read file"))
    r.readAsText(file)
})

const parseVariantsFromText = (text) => text.split(/\r?\n/).map((l) => l.trim()).filter((l) => l && !l.startsWith("#"))

const parseVariantsFromVCF = (content) => {
    const variants = []
    for (const line of content.split(/\r?\n/)) {
        if (!line.trim() || line.startsWith("#")) continue
        const cols = line.split("\t")
        if (cols.length < 5) continue
        const chrom = cols[0].replace(/^chr/i, "")
        const pos = cols[1]
        const ref = cols[3]
        const altField = cols[4]
        if (ref === "." || altField === "." || !ref || !altField) continue
        if (ref.includes("<") || altField.includes("<")) continue
        const alt = altField.split(",")[0].trim()
        if (!alt) continue
        variants.push(`${chrom}-${pos}-${ref}-${alt}`)
    }
    return variants
}

const resolveVariantListFromInputs = async () => {
    const fileInput = document.getElementById("batch-vcf-file")
    if (fileInput && fileInput.files && fileInput.files.length > 0) {
        const file = fileInput.files[0]
        const name = file.name.toLowerCase()
        if (name.endsWith(".vcf.gz")) {
            throw new Error("Please upload a plain-text .vcf file (gzip is not supported in the browser).")
        }
        const text = await readFileAsText(file)
        if (name.endsWith(".vcf")) {
            return parseVariantsFromVCF(text)
        }
        return parseVariantsFromText(text)
    }
    const batchText = $("#batch-variants-textarea").val().trim()
    if (batchText) {
        const lines = parseVariantsFromText(batchText)
        if (lines.length) return lines
    }
    const single = $("#search-box").val().trim()
    return single ? [single] : []
}

const initResultTablePopups = () => {
    $(".question, .exclamation, .score-table").popup({
        "on": "click",
        "lastResort": "bottom left",
    })
    $(".question").css({"color": "#666666", "cursor": "pointer"})
    $(".exclamation").css({"cursor": "pointer"})
    $(".score-table").css({"color": "#000080", "cursor": "pointer"})
    $(".ui.modal").modal()
}

const refreshBatchNavControls = () => {
    const n = batchVariantResults.length
    $("#batch-prev-btn").prop("disabled", currentBatchIndex <= 0)
    $("#batch-next-btn").prop("disabled", currentBatchIndex >= n - 1 || n === 0)
}

const displayBatchVariantAtIndex = async (idx) => {
    if (!batchVariantResults.length || idx < 0 || idx >= batchVariantResults.length) return
    if (!lastBatchFormOptions) return

    currentBatchIndex = idx
    const o = lastBatchFormOptions
    const entry = batchVariantResults[idx]

    $("#batch-counter").text(`Variant ${idx + 1} of ${batchVariantResults.length}`)
    $("#batch-variant-select").val(String(idx))
    refreshBatchNavControls()

    $("#error-box").html("")

    if (entry.normalizeError) {
        showError(entry.normalizeError)
        $("#spliceai-table, #pangolin-table, #transcript-button-table, #other-predictors-table").hide()
        $("#response-box").show()
        return
    }

    const variant = entry.rawVariant
    const normalizedVariant = entry.normalizedVariant
    const variantConsequence = entry.consequence

    lastSpliceaiResponseJson = entry.spliceaiJson
    lastPangolinResponseJson = entry.pangolinJson
    lastGenomeVersion = o.genomeVersion

    $(".ref-score-column, .alt-score-column").toggle(o.showRefAltColumns == "1")

    if (!entry.spliceaiError && entry.spliceaiJson) {
        renderSplicingResultsFromApiJson(entry.spliceaiJson, normalizedVariant, variantConsequence, "SpliceAI", variant, o.genomeVersion, o.basicOrComprehensive, o.maxDistance, o.mask, o.showRefAltColumns)
        $("#spliceai-table, #transcript-button-table").show()
    } else {
        if (entry.spliceaiError) showError(entry.spliceaiError)
        $("#spliceai-table, #transcript-button-table").hide()
    }

    if (!entry.pangolinError && entry.pangolinJson) {
        renderSplicingResultsFromApiJson(entry.pangolinJson, normalizedVariant, variantConsequence, "Pangolin", variant, o.genomeVersion, o.basicOrComprehensive, o.maxDistance, o.mask, o.showRefAltColumns)
        $("#pangolin-table").show()
    } else {
        if (entry.pangolinError) showError(entry.pangolinError)
        $("#pangolin-table").hide()
    }

    await generateOtherPredictorsTable(normalizedVariant, variantConsequence, variant, o.genomeVersion)
    $("#other-predictors-table").show()

    $("#response-box").show()
    initResultTablePopups()

    if (lastSpliceaiResponseJson != null || lastPangolinResponseJson != null) {
        $("#igv-table").show()
        $("#igv-div").hide()
        $("#show-igv-button").text("Show")
    }

    window.location.hash = "#" + $.param({
        variant: variant,
        hg: o.genomeVersion,
        bc: o.basicOrComprehensive,
        distance: o.maxDistance,
        mask: o.mask,
        ra: o.showRefAltColumns,
    })
}

const runBatchVariantSubmit = async (variants, formOptions) => {
    lastBatchFormOptions = formOptions
    batchVariantResults = []
    $("#batch-nav").show()
    $("#batch-progress").text("")

    for (let i = 0; i < variants.length; i++) {
        $("#batch-progress").text(`Loading variant ${i + 1} of ${variants.length}…`)
        const rawVariant = variants[i]
        const entry = {
            rawVariant,
            normalizedVariant: null,
            consequence: null,
            spliceaiJson: null,
            pangolinJson: null,
            spliceaiError: null,
            pangolinError: null,
            normalizeError: null,
        }
        try {
            const norm = await normalizeVariant(rawVariant, formOptions.genomeVersion)
            entry.normalizedVariant = norm.variant
            entry.consequence = norm.consequence
        } catch (e) {
            entry.normalizeError = e.message
            batchVariantResults.push(entry)
            continue
        }

        const [spliceaiSettled, pangolinSettled] = await Promise.allSettled([
            fetchSplicingToolJson(entry.normalizedVariant, entry.consequence, "SpliceAI", rawVariant, formOptions.genomeVersion, formOptions.basicOrComprehensive, formOptions.maxDistance, formOptions.mask),
            fetchSplicingToolJson(entry.normalizedVariant, entry.consequence, "Pangolin", rawVariant, formOptions.genomeVersion, formOptions.basicOrComprehensive, formOptions.maxDistance, formOptions.mask),
        ])

        if (spliceaiSettled.status === "fulfilled") {
            entry.spliceaiJson = spliceaiSettled.value
        } else {
            entry.spliceaiError = spliceaiSettled.reason.message || String(spliceaiSettled.reason)
        }
        if (pangolinSettled.status === "fulfilled") {
            entry.pangolinJson = pangolinSettled.value
        } else {
            entry.pangolinError = pangolinSettled.reason.message || String(pangolinSettled.reason)
        }

        batchVariantResults.push(entry)
    }

    $("#batch-progress").text("")
    const $sel = $("#batch-variant-select").empty()
    batchVariantResults.forEach((e, i) => {
        const short = e.rawVariant.length > 72 ? `${e.rawVariant.slice(0, 69)}…` : e.rawVariant
        const label = e.normalizeError ? `${short} (failed)` : short
        $sel.append($("<option/>").attr("value", i).text(`${i + 1}. ${label}`))
    })

    await displayBatchVariantAtIndex(0)
}

const runSingleVariantSubmit = async (variant, formOptions) => {
    const genomeVersion = formOptions.genomeVersion
    const basicOrComprehensive = formOptions.basicOrComprehensive
    const maxDistance = formOptions.maxDistance
    const mask = formOptions.mask
    const showRefAltColumns = formOptions.showRefAltColumns

    try {
        const normalizedVariantAndConsequence = await normalizeVariant(variant, genomeVersion)
        const normalizedVariant = normalizedVariantAndConsequence.variant
        const variantConsequence = normalizedVariantAndConsequence.consequence

        let [spliceaiResponseJson, pangolinResponseJson, otherPredictorsResponseJson] = await Promise.all([
            generateSplicingResultsTable(normalizedVariant, variantConsequence, "SpliceAI", variant, genomeVersion, basicOrComprehensive, maxDistance, mask, showRefAltColumns),
            generateSplicingResultsTable(normalizedVariant, variantConsequence, "Pangolin", variant, genomeVersion, basicOrComprehensive, maxDistance, mask, showRefAltColumns),
            generateOtherPredictorsTable(normalizedVariant, variantConsequence, variant, genomeVersion),
        ].map(p => p.catch(e => ({'error': e.message}))))

        $("#response-box").show()

        initResultTablePopups()

        if (spliceaiResponseJson.error) {
            showError(`${spliceaiResponseJson.error}`)
            spliceaiResponseJson = null
            $("#spliceai-table, #transcript-button-table").hide()
        } else {
            $("#spliceai-table, #transcript-button-table").show()
        }

        if (pangolinResponseJson.error) {
            showError(`${pangolinResponseJson.error}`)
            pangolinResponseJson = null
            $("#pangolin-table").hide()
        } else {
            $("#pangolin-table").show()
        }

        window.location.hash = "#" + $.param({
            variant: variant,
            hg: genomeVersion,
            bc: basicOrComprehensive,
            distance: maxDistance,
            mask: mask,
            ra: showRefAltColumns,
        })

        lastGenomeVersion = genomeVersion
        lastSpliceaiResponseJson = spliceaiResponseJson
        lastPangolinResponseJson = pangolinResponseJson

        if (spliceaiResponseJson != null || pangolinResponseJson != null) {
            $("#igv-table").show()
            $("#igv-div").hide()
            $("#show-igv-button").text("Show")
        }
    } catch(e) {
        console.error(e)
        showError(e.message)
    }
}

const updateIgvBrowser = async (spliceaiResponseJson, pangolinResponseJson, genomeVersion) => {
    if ((spliceaiResponseJson && spliceaiResponseJson.allNonZeroScores) || (
        pangolinResponseJson && pangolinResponseJson.allNonZeroScores)) {
        const igvConfig = generateIgvConfig(spliceaiResponseJson, pangolinResponseJson, genomeVersion)
        if (!window.igvBrowser) {
            // create igv.js browser object if it hasn't been created yet
            window.igvBrowser = await igv.createBrowser(document.getElementById("igv-viewport"), igvConfig)
        } else {
            console.log("Updating igv config to", igvConfig)
            await window.igvBrowser.loadSessionObject(igvConfig)
        }
    }
}


const getFormOptions = () => ({
    genomeVersion: $("input[name='hg']:checked").val().trim(),
    basicOrComprehensive: $("input[name='gencode-gene-set']:checked").val().trim(),
    maxDistance: $("#max-distance-input").val().trim(),
    mask: $(`input[name='mask']`).prop("checked") ? "1" : "0",
    showRefAltColumns: $(`input[name='show-ref-alt']`).prop("checked") ? "1" : "0",
})

const handleSubmit = async () => {
    const formOptions = getFormOptions()
    const genomeVersion = formOptions.genomeVersion

    updateVisualizationCheckboxes(false, genomeVersion)

    $("#submit-button").addClass("loading disabled")
    $("#response-box, #spliceai-table, #pangolin-table, #transcript-button-table, #error-box").hide()
    $("#error-box").html("")

    let variants = []
    try {
        variants = await resolveVariantListFromInputs()
    } catch (e) {
        showError(e.message)
        $("#submit-button").removeClass("loading disabled")
        return
    }

    if (!variants.length) {
        showError("No variants to analyze. Enter a variant in the top field, paste one or more lines in the batch box, or upload a plain-text .vcf file.")
        $("#submit-button").removeClass("loading disabled")
        return
    }

    if (variants.length > BATCH_VARIANT_MAX) {
        showError(`Please enter at most ${BATCH_VARIANT_MAX} variants per batch.`)
        $("#submit-button").removeClass("loading disabled")
        return
    }

    $(".ref-score-column, .alt-score-column").toggle(formOptions.showRefAltColumns == "1")

    if (variants.length === 1) {
        $("#batch-nav").hide()
        batchVariantResults = []
        lastBatchFormOptions = null
        $("#batch-progress").text("")
        await runSingleVariantSubmit(variants[0], formOptions)
        $("#submit-button").removeClass("loading disabled")
        return
    }

    await runBatchVariantSubmit(variants, formOptions)
    $("#submit-button").removeClass("loading disabled")
}

const applyUrlSettingsToFormElements = () => {
    // get optional settings from the url and update form settings
    let hgFromUrl = $.urlParam('hg')
    if (hgFromUrl == "19") {
        hgFromUrl = "37"
    }
    if (hgFromUrl && (hgFromUrl == "38" || hgFromUrl == "37")) {
        $(`input[name='hg'][value='${hgFromUrl}']`).prop("checked", true)
    }

    const bcFromUrl = $.urlParam('bc')
    if (bcFromUrl && (bcFromUrl == "basic" || bcFromUrl == "comprehensive")) {
        $(`input[name='gencode-gene-set'][value='${bcFromUrl}']`).prop("checked", true)
    }

    const maxDistanceFromUrl = $.urlParam('distance')
    if (maxDistanceFromUrl) {
        $("#max-distance-input").val(maxDistanceFromUrl)
    }

    const maskFromUrl = $.urlParam('mask')
    if (maskFromUrl) {
        $(`input[name='mask']`).prop("checked", maskFromUrl == "1")
    }

    const showRefAltColumnsFromUrl = $.urlParam('ra')
    if (showRefAltColumnsFromUrl) {
        $(`input[name='show-ref-alt']`).prop("checked", showRefAltColumnsFromUrl == "1")
    }

    //update the variant input box last and trigger a search
    const variantFromUrl = $.urlParam('variant')
    if (variantFromUrl) {
        $("#search-box").val(variantFromUrl)
        $("#submit-button").click()
    }
}

/*
const toggleRefAltScoreColumns = () => {
    const currentText = $("#toggle-ref-alt-scores-button").html()

    if (currentText.match(/Show/i)) {
        $(".ref-alt-score-column").show()
        $("#toggle-ref-alt-scores-button").removeClass('primary')
        $("#toggle-ref-alt-scores-button").html('Hide REF & ALT Scores')
    } else {
        $(".ref-alt-score-column").hide()
        $("#toggle-ref-alt-scores-button").addClass('primary')
        $("#toggle-ref-alt-scores-button").html('Show REF & ALT Scores')
    }
}
*/

// define function for parsing url parameters like ?variant=chr1-1234567-T-G
$.urlParam = (name) => {
    const results = new RegExp('[\?&#]?' + name + '=([^&#]*)').exec(window.location.hash);
    return (results !== null) ? decodeURIComponent(results[1]) || 0 : false;
}

$(document).ready(() => {
    $("#gencode-version").text(GENCODE_VERSION)
    // init UI elements
    $(".ui.checkbox").checkbox()
    $(".question, .exclamation, .score-table").popup({"on": "click"})
    $(".question, .exclamation, .score-table").css("cursor", "pointer")
    $("#search-box").focus()
    $("#response-box").hide()

    // init event handlers
    $("#submit-button").click(handleSubmit)

    $("#batch-prev-btn").click(() => { void displayBatchVariantAtIndex(currentBatchIndex - 1) })
    $("#batch-next-btn").click(() => { void displayBatchVariantAtIndex(currentBatchIndex + 1) })
    $("#batch-variant-select").on("change", function() {
        void displayBatchVariantAtIndex(parseInt($(this).val(), 10))
    })

    $("#search-box, #max-distance-input").keydown((event) => {
        if ((event.keyCode || event.which) == 13) {
            $("#submit-button").click()
        }
    })

    $("#show-igv-button").click(async () => {
        $("#show-igv-button").text("Update")
        $("#igv-div").show()
        await updateIgvBrowser(lastSpliceaiResponseJson, lastPangolinResponseJson, lastGenomeVersion)
    })

    /*
    $("#show-other-predictors-button").click(() => {
        $(".other-predictors-row").show()
        $("#show-other-predictors-row").hide()
    })
    */

    updateVisualizationCheckboxes(true, "38")

    applyUrlSettingsToFormElements()
})