#imports: DFU, others....
from installed_clients.DataFileUtilClient import DataFileUtil
from copy import deepcopy

def DownloadMotifSet(refList,callback):
    MotifSetDict = {}
    #init DFU
    dfu = DataFileUtil(callback)
    for ref in refList:
        get_objects_params = {'object_refs' : [ref]}
        #get_ss_params = {'object_refs' : [params['SS_ref']]}
        MotifSet = dfu.get_objects(get_ss_params)['data'][0]['data']
        MotifSetDict[ref] = deepcopy(MotifSet)
    return MotifSetDict
