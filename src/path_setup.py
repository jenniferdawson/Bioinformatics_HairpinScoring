

BASE_PATH = ""  # ie: /home/user/HairpinScoring/
DSSP_PATH = ""  # Optional (download folder for third party DSSP database)

if BASE_PATH == "":
    print("ERROR: Set BASE_PATH in ./src/path_setup.py")
    exit()

DB_PATH = BASE_PATH+"/database/"
DERIVED_PATH = BASE_PATH+"/derived/"
PLOT_PATH = BASE_PATH+"/plots/"