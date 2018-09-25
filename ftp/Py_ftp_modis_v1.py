# ftp-ex.py
# http://www.blog.pythonlibrary.org/2012/07/19/python-101-downloading-a-file-with-ftplib/
# http://www.pythonforbeginners.com/code-snippets-source-code/how-to-use-ftp-in-python
# https://docs.python.org/2/library/ftplib.html
# https://docs.python.org/2/tutorial/errors.html

import sys
sys.path.insert(0, 'C:\\Users\\carlos\\Dropbox\\newPython\\utils')
sys.path.insert(0, 'C:\\Users\\carlos.silveira\\Dropbox\\newPython\\utils')

import os
from ftplib import FTP
import os.path
import tarfile
import gc
from utils import printlog, sleepTime, mapDict
from time import strftime, localtime
from sys import argv
import threading
import queue
import time
import msvcrt

usage = """\
Usage: %s [OPTIONS]
        -fo     folder to output data
        -ds     date start in year|doy
        -de     date end in year|doy
        -ty     type of data [MOD09GA, MYD09GA ]
                See more types in ftp://ladsweb.nascom.nasa.gov/
""" % argv[0]

threadList = ["Thread-1", "Thread-2", "Thread-3", "Thread-4", "Thread-5"]
queueLock = threading.Lock()
workQueue = queue.Queue(10)
threads = []
exitFlag = 0

# date is year + doy
def ftpCopy(outputDir, date, typeData):
    date = str(date)
    year = date[0:4]
    doy = date[4:7]
    try:
        ftp = FTP('ladsweb.nascom.nasa.gov', 'anonymous')
    except:
        #print( 'Error! Don''t have Connection to Internet.')
        return ([], -1)
    try:
        ftp.login()
        print('Logged in ftp server')
        #print(ftp.getwelcome())
    except:
        if(ftp.getwelcome()== None):
            print('Error! Verify user name and password')
            return ([], -2)
        #print('Finally in login')
    files = ['.h12v10.006', '.h13v10.006']
    try:
        ftp.cwd("allData")
        ftp.cwd("6")
        ftp.cwd(typeData)
        ftp.cwd(year)
        ftp.cwd(doy)
    except:
        print('Error in ftp folder: %s %s'%(year, doy))
        # Test if don't get data in folder
        if(typeData in ['MOD09GA', 'MYD09GA']):
            return ([], -3)
        else:
            return ([], 0)
    # ftp.cwd("subFolder")
    # or ftp.cwd("folderOne\\subFolder")
    #ftp.retrlines("LIST")
    listing = []
    lstFiles = []
    ftp.retrlines("LIST", listing.append)
    poslist = -1
    filename = []
    cont_file = 0
    #print(listing)
    for n_file in range(len(files)):
        try:
            while True:
                poslist = poslist + 1
                #print('.... File: ' + str(poslist))
                words = listing[poslist].split(None, 8)
                filename = words[-1].lstrip()
                # Filter to year data
                #print(filename)
                if(files[0] in str(filename) or  files[1] in str(filename)):
                    cont_file = 1
                    break
        except:
            if (poslist > len(listing)):
                print('Finish copy files')
                cont_file = 0
                break
        if(cont_file == 1):
            cont_file = 0
            # download the file
            subfolder = year + '\\' + doy + '\\'
            #print(subfolder)
            if not(os.path.isdir(outputDir + subfolder)):
                os.makedirs(outputDir + subfolder)
            try:
                if(os.path.isfile(outputDir + subfolder + filename) == True):
                    print( 'File yet copy: ' + outputDir + subfolder + filename)
                else:
                    print( 'File to copy: ' +  outputDir + subfolder + filename)
                    lstFiles.append(filename)
            except:
                print('.......Error in names (wait):')
                print(outputDir)
                print(subfolder)
                print(filename)
                msvcrt.getch()
        else:
            break

    ftp.quit()
    return (lstFiles, 0)

class myThread (threading.Thread):
    def __init__(self, threadID, name, q, outputDir, typeData):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.q = q
        self.outputDir = outputDir
        self.typeData = typeData
    def run(self):
        printlog (True, "Starting " + self.name)
        #print( "q: " + str(self.q) + " typeData: " + self.typeData + " donwloadDir: " + self.outputDir)
        process_data(self.name, self.q, self.outputDir, self.typeData)
        printlog (True, "Exiting " + self.name)

def process_data(threadName, q, outputDir, typeData):
    global queueLock
    global workQueue
    global exitFlag
    #print( "Enter in process_data. exitFlag: " + str(exitFlag))

    while not exitFlag:
        #print( "exitFlag is False")
        queueLock.acquire()
        #print( "\nqueueLock acquire")
        if not workQueue.empty():
            #print( "Queue is not empty")
            filename = q.get()
            queueLock.release()
            #print( "\nqueueLock release")
            print ("\n%s processing %s" % (threadName, filename))
            error = -5
            while(error != 0):
                error = copyFTPfile(filename, outputDir, typeData)
                if(error != 0):
                    print('... error %d wait thread %s'%(error, str(threadName)))
                    time.sleep(10)
        else:
            #print( "Queue is empty")
            try:
                queueLock.release()
                #print( "\nqueueLock release")
            except:
                print( 'Fail to release queueLock !')
        time.sleep(1)
    print( "exitFlag is True")

def copyFTPfile(inputFile, outputDir, typeData):
    try:
        ftp = FTP('ladsweb.nascom.nasa.gov', 'anonymous')
    except:
        print( 'Error! Don''t have Connection to Internet.')
        return(-1)
    try:
        ftp.login()
        #print('Logged in ftp server')
        #print(ftp.getwelcome())
    except:
        if(ftp.getwelcome()== None):
            print('Error! Verify user name and password')
            return(-2)
        #print('Finally in login')
    if(typeData in ['MOD09GA', 'MYD09GA', 'MOD13A2', 'MOD13Q1']):
        year  = inputFile[9:13]
        doy = inputFile[13:16]
    else:
        year  = inputFile[10:14]
        doy = inputFile[14:17]
    try:
        ftp.cwd("allData")
        ftp.cwd("6")
        ftp.cwd(typeData)
        ftp.cwd(year)
        ftp.cwd(doy)
    except:
        print('Error in ftp folder: %s%s'%(year, doy))
        return (-3)
    subfolder = year + '\\' + doy + '\\'
    if(typeData in ['MOD09GA', 'MYD09GA']):
        outputFile = outputDir + subfolder + inputFile
    else:
        outputFile = outputDir + inputFile
    print('. Save file in: %s'%(outputFile))
    #printlog(True, "Create FTP data file")
    error = 0
    while(error < 5):
        try:
            lf = open(outputFile, "wb")
            ftp.retrbinary("RETR " + inputFile, lf.write, 8*1024)
            lf.close()
            return (0)
        except:
            error = error + 1
        if(error >= 5):
            print("Problem to copy file: %s"%(outputFile))
            return (-4)    

class main:

    global exitFlag
    global threads
    global workQueue
    threadID = 1

    argDict = mapDict(argv, usage)

    if("-fo" in argDict and "-ds" in argDict and "-de" in argDict and "-ty" in argDict):
        outputDir = argDict["-fo"]
        dateStart   = argDict["-ds"]
        dateFinish  = argDict["-de"]
        typeData = str(argDict["-ty"])
    else:
        exit(usage)

    outputDir = outputDir + '\\'     
    gc.collect()
    # Create file list to download
    lst_all = []
    date = int(dateStart)
    while(date < int(dateFinish)):
        print('Verify files do date: %i'%date)
        process = False
        while(process == False):
            [lst, error] = ftpCopy(outputDir, date, typeData)
            print(lst)
            print(error)
            if(error == 0):
                process = True
            else:
                print('...')
                time.sleep(10)
        for id in range(len(lst)):
            lst_all.append(lst[id])
        #print('lst_all: ' + str(len(lst_all)))
        #print(lst_all)
        # Increment date [year|doy]
        year = int(str(date)[0:4])
        doy = int(str(date)[4:7])
        if doy < 367:
            if(typeData in ['MCD15A2H', 'MOD15A2H']):
                doy = doy + 8
            else:
                if(typeData in ['MOD13Q1', 'MOD13A2']):
                    doy = doy + 16
                else:
                    doy = doy + 1
        else:
            year = year + 1
            doy = 1
        date = int(str('%04d'%year) + str('%03d'%doy))

    # download file list
    print('-------- lst_all ----------')
    print(lst_all)
    if(lst_all == []):    
        print( 'Wait new files')
        sys.exit()

    print( "----- Create new threads")
    for tName in threadList:
        thread = myThread(threadID, tName, workQueue, outputDir, typeData)
        thread.start()
        threads.append(thread)
        threadID += 1

    print( "----- Fill the queue and wait case many files")
    
    x = 0
    while(x < len(lst_all)):
        queueLock.acquire()
        if not workQueue.full():
            filename =  lst_all[x]
            print( 'x: ' + str(x) + '\tfilename: ' + filename)
            workQueue.put(filename)
            x = x + 1
            
        else:
            time.sleep(30)
        queueLock.release()

    print( "----- Wait for queue to empty")
    while not workQueue.empty():
        #print('Files in Queue: ' + str(workQueue.qsize()))
        #time.sleep(5)
        pass
    print( '----- Queue is empty')

    # Notify threads it's time to exit
    exitFlag = 1

# Wait for all threads to complete
    for t in threads:
        t.join()
    print( "Exiting Main Thread")
    print( 'Waiting new files')
    
