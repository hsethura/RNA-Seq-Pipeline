#!/usr/bin/env python

from __future__ import print_function
import httplib2
import os

from apiclient import discovery
from oauth2client import client
from oauth2client import tools
from oauth2client.file import Storage

try:
    import argparse
    parser = argparse.ArgumentParser(parents=[tools.argparser])
    parser.add_argument('--spid', '-s', type = str, required = True, help = 'Google spreadsheet id')
    parser.add_argument('--tabname', '-t', type = str, required = False, default = "Sample Information", help = 'Tab/Sheet name in the file')
    parser.add_argument('--proid', '-p', type = str, required = True, help = 'Project Id to generate key file')
    parser.add_argument('--Key_dir', type = str, required = False, default = 'none', help = 'Key_dir to store downloaded key files')
    flags = parser.parse_args()
except ImportError:
    flags = None

# If modifying these scopes, delete your previously saved credentials
# at ~/.credentials/sheets.googleapis.com-python-quickstart.json
SCOPES = 'https://www.googleapis.com/auth/spreadsheets.readonly'
CLIENT_SECRET_FILE = 'client_secret.json'
APPLICATION_NAME = 'Google Sheets API Python Quickstart'


def get_credentials():
    """Gets valid user credentials from storage.

    If nothing has been stored, or if the stored credentials are invalid,
    the OAuth2 flow is completed to obtain the new credentials.

    Returns:
        Credentials, the obtained credential.
    """
    home_dir = os.path.expanduser('~')
    credential_dir = os.path.join(home_dir, '.credentials')
    if not os.path.exists(credential_dir):
        os.makedirs(credential_dir)
    credential_path = os.path.join(credential_dir,
                                   'sheets.googleapis.com-python-quickstart.json')

    store = Storage(credential_path)
    credentials = store.get()
    if not credentials or credentials.invalid:
        flow = client.flow_from_clientsecrets(CLIENT_SECRET_FILE, SCOPES)
        flow.user_agent = APPLICATION_NAME
        if flags:
            credentials = tools.run_flow(flow, store, flags)
        else: # Needed only for compatibility with Python 2.6
            credentials = tools.run(flow, store)
        print('Storing credentials to ' + credential_path)
    return credentials

def main():
    """Shows basic usage of the Sheets API.

    Creates a Sheets API service object and prints the names and majors of
    students in a sample spreadsheet:
    https://docs.google.com/spreadsheets/d/1BxiMVs0XRA5nFMdKvBdBZjgmUUqptlbs74OgvE2upms/edit
    """
    credentials = get_credentials()
    http = credentials.authorize(httplib2.Http())
    discoveryUrl = ('https://sheets.googleapis.com/$discovery/rest?'
                    'version=v4')
    service = discovery.build('sheets', 'v4', http=http,
                              discoveryServiceUrl=discoveryUrl)
    spreadsheetId = flags.spid
    tabname =flags.tabname
    #spreadsheetId = '1R3ZjJjevoqoNKagR-7xEShsYQ1O-DQILfTOP1rANItc'
    rangeName = tabname
    result = service.spreadsheets().values().get(
        spreadsheetId=spreadsheetId, range=rangeName).execute()
    values = result.get('values', [])
    #values = result.get('values')

    #outfile = 'outfile.tsv'
    #outfile = flags.outfile
    proid = flags.proid

    KeyDir = ''
    if flags.Key_dir == 'none':
        KeyDir = "/broad/IDP-Dx_storage/MOC/Key_files"
    else:
        KeyDir = flags.Key_dir

    outfile = KeyDir + "/" + proid  + "_key.txt"
    outf = open(outfile, 'w')

    seen_header = False
    if not values:
        print('No data found.')
    else:
        #print('Name, Major:')
        header_len = 0
        for lrow in values:
            if not lrow[0].startswith("###"):
                if not seen_header:
                    seen_header = True
                    header_len = len(lrow)
                else:
                    lrow_len = len(lrow)
                    len_diff = header_len - lrow_len
                    #print("len_diff: " + str(len_diff))
                    tab_lst = list("\t" * len_diff)
                    lrow = lrow + tab_lst
                    
  
            count = 0
            #print(': '.join(lrow))
            lrow_len = len(lrow)
            for lvalue1 in lrow:
                lvalue = lvalue1.rstrip()
                if count == 0:
                    outf.write(lvalue.encode('utf-8'))
                    count = count + 1
                else:
                    outf.write('\t'.encode('utf-8'))
                    outf.write(lvalue.encode('utf-8'))
            outf.write('\n')

	print(outfile)
if __name__ == '__main__':
    main()

