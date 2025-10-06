import numpy as np
# import pandas as pd
from pyvo.dal import TAPService
from pyvo.dal.adhoc import DatalinkResults
# from astroquery.alma import Alma
from time import sleep

loc = {'US': "almascience.nrao.edu", 'EU':"almascience.eso.org"}

def set_service(location):
	return TAPService(f"https://{loc[location]}/tap")

def server_choice(n):
    if((n%2)==0):
        return 'US'
    else:
        return 'EU'

def server(n):
    return set_service(server_choice(n))

def query_proposal_id(service, proposal_id):
	""" QUeries for an ALMA project code
		service 	pyvo TAPService instance
		project_id 	ALMA proposal id in the form "2012.1.00123.S" or in the form "ADS/JAO.ALMA#2012.1.00123.S"

		return 		pandas Dataframe"""

	query = f"""
			SELECT *
			FROM ivoa.obscore
			WHERE obs_publisher_did like '%{proposal_id}%'
			"""
	return service.search(query).to_table().to_pandas()


def query_member_ids(service, proposal_id):
  	return np.unique(query_proposal_id(service, proposal_id)['member_ous_uid'])

def query_datalink_core(location, member_uid):
    query_link = f"https://{loc[location]}/datalink/sync?ID={member_uid}"
    print(query_link)
    datalink = DatalinkResults.from_result_url(query_link)
    return datalink.to_table().to_pandas()


def query_datalink(location, member_uid, retry_time=5, sleep_time=10):
    connection_time = 0
    temp = None

    while connection_time < retry_time:
        try:
            temp = query_datalink_core(location, member_uid)
            if connection_time > 0:
                print('Re-connection successful')
            break
        except:
            print(f'{member_uid} connection LOST')
            connection_time += 1
            sleep(sleep_time)

    if not temp:
        print("Connection failed")

    return temp

def mous_query(pid):
    done = False
    n = 0
    mous = ''
    while not done:
        serv = server(n)

        try:
            mous = query_member_ids(serv, pid)
            sleep(5)
            done = True
        except: # DALFormatError
            if (n>=20):
                mous = ''
                n = 100
                done = True
            print("CONNECT LOST, ... RECONNECTIONG")
            sleep(10)
            n += 1
    return mous, n

def size_query(uid):
    serv_list = ['US', 'EU']
    done = False
    n = 0

    sq= ''

    while not done:
        serv = serv_list[n%2]
        try:
            sq = query_datalink(serv, uid)
            done = True
            sleep(10)
        except:
            if(n>=20):
                sq = ''
                n = 100
                done = True
            print('CONNECT LOST, ... RECONNECTIONG')
            n += 1
            sleep(10)
    return sq, n


