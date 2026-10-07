from kolach.metabolism.dbcan import download_dbcan, DATABASE_FILES

rule download_dbcan_references:
    output:
        manifest=f"{DBCAN_DIR}/kolach_dbcan_database.json",
        references=[f"{DBCAN_DIR}/{name}" for name in DATABASE_FILES]
    run:
        download_dbcan(DBCAN_DIR)
