/*
 * HTTP entry point for the score-lookup Cloud Run function.
 *
 * GET /?hg=38&chrom=8&pos=140300616&ref=T&alt=G
 *   200 {"variant": {...}, "scores": {"alphagenome_avi": {"rawScore": 1.59, "phredScore": 30.48}, ...}, "errors": []}
 *   400 {"error": "..."} when a parameter is missing or malformed
 *
 * "scores" holds only the scores found for the variant. "errors" lists any table whose query failed,
 * in which case the response isn't cached.
 */
import functions from '@google-cloud/functions-framework'

import { createTables, lookupScores, parseVariantParams } from './lookup.js'

// Created once per instance, so each table's index is downloaded and parsed on its first query and
// then reused by every later request the instance serves.
const tables = createTables()

// An hour covers re-running the same search, and the request URL doesn't change when a deploy
// switches to new tables or parsing, so a longer cache would keep serving the old answer that much longer.
const CACHE_CONTROL_SUCCESS = 'public, max-age=3600'

export const handleRequest = async (req, res, tablesToQuery = tables) => {
    res.set('Access-Control-Allow-Origin', '*')
    if (req.method === 'OPTIONS') {
        res.set('Access-Control-Allow-Methods', 'GET')
        res.set('Access-Control-Max-Age', '3600')
        res.status(204).send('')
        return
    }
    if (req.method !== 'GET') {
        res.set('Allow', 'GET, OPTIONS')
        res.status(405).json({ error: 'Only GET requests are supported' })
        return
    }

    const { variant, error } = parseVariantParams(req.query)
    if (error) {
        res.set('Cache-Control', 'no-store')
        res.status(400).json({ error })
        return
    }

    const { scores, errors } = await lookupScores(tablesToQuery, variant)
    for (const tableError of errors) {
        console.error(`Query of ${tableError.table} failed for ${JSON.stringify(variant)}: ${tableError.message}`)
    }
    res.set('Cache-Control', errors.length ? 'no-store' : CACHE_CONTROL_SUCCESS)
    res.status(200).json({ variant, scores, errors })
}

functions.http('scoreLookup', (req, res) => handleRequest(req, res))
