import sqlite3


def load_database(card, name):
    """
    Build a database to contain
    """

    db_fp = "{}_training.db".format(name)
    if os.path.exists(db_fp):
        print("{} db already exists, loading.".format(db_fp))
    else:
        create_database(card, db_fp)

    conn = sqlite3.connect(db_fp)

    return conn


def create_database(card, db_fp):
    """
    Create the database from the parsed CARD data
    Primary key is ARO

    Table of ARO with family, nucleotide sequence, and protein sequence

    This table has a one-to-many relationship with table 2:
    Table of overlaps within family
    ARO, ARO2, overlaps of readlength
    """

    conn = sqlite3.connect(db_fp)

    c = conn.cursor()

    c.execute('''CREATE TABLE sequences
                    (aro text, family text, nt text, aa text)''')

    c.execute('''CREATE TABLE overlaps
                    (aro text, aro2 text, size int, overlaps int)''')



