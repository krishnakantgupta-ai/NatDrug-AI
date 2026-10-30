# Before a public GitHub push

- [ ] Copy only the supplied `research_v2/` folder into the existing `NatDrug-AI` repository. Keep the existing `app.py`, README, and any project LICENSE intact.
- [ ] Do not publish the **local-only data archive** or confidential ACS reviewer comments, manuscripts, email, institutional credentials, or unpublished private source records.
- [ ] Confirm your repository's current software license (none is supplied or changed here).
- [ ] Give ChEMBL attribution and retain compatible share-alike treatment for redistributed ChEMBL-derived datasets.
- [ ] Confirm IMPPAT public redistribution rights before uploading bulk exports or full structure-level derived screening CSVs.
- [ ] Run `python -m unittest discover -s research_v2/tests -v` from the repository root after installing pandas.
- [ ] Report measured CV and in-silico predictions distinctly. Do not claim independent natural-product accuracy until original records are verified.
- [ ] If you want a permanent DOI, create a versioned GitHub Release and archive that version with Zenodo after reviewing included assets.
