
import pyodbc

server = 'UEORS-DB\\'
database = 'UEORS_MAIN_DB'
# cnxn = pyodbc.connect('DRIVER={ODBC Driver 17 for SQL Server};SERVER=' + server + ';DATABASE=' + database + ';Trusted_Connection=yes;')
username = 'ueors_user'
password = 'ueor101'
cnxn = pyodbc.connect('DRIVER={ODBC Driver 17 for SQL Server};SERVER='+server+';DATABASE='+database+
                      ';UID='+username+';PWD='+ password)
cursor = cnxn.cursor()

try:
    result = cursor.execute('SELECT * FROM companies')
    print(result.fetchall())
except Exception as e:
    print(e)
finally:
    cnxn.close()
