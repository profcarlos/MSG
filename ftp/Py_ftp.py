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

usage = """\
Usage: %s [OPTIONS]
        -fo     folder to output data
        -ty     type input data (in FTP server: GBR, HRIT, HDF)
""" % argv[0]

threadList = ["Thread-1", "Thread-2", "Thread-3"] #, "Thread-4", "Thread-5"]
queueLock = threading.Lock()
workQueue = queue.Queue(10)
threads = []
exitFlag = 0

def ftpCopy(downloadDir, typeData):
    try:
        ftp = FTP("ftp.lapig.iesa.ufg.br", "lapig-ftp", "lapig123ftp")
    except:
        printlog(True, 'Error! Don''t have Connection to Internet.')
        return(-1)
    try:
        ftp.login()
        printlog(True, 'Logged in ftp server')
        printlog(ftp.getwelcome())
    except:
        if(ftp.getwelcome()!= None):
            printlog(True, 'Logged in ftp server')
            printlog(True, ftp.getwelcome())
        else:
            printlog(True, 'Error! Verify user name and password')
            return(-2)
        printlog(False, 'Finally in login')

    # To cloud mask filtes
    if(typeData == 'GRB'):
        filter1 = 'MSG3-SEVI-MSGCLMK-0100-0100'
        filter2 = '.grb'
    elif(typeData == 'HRIT'):
        filter1  = 'MSG3-SEVI-MSG15-0100-NA-201'
        filter2  = '.tar'
    elif(typeData == 'HDF'):
        filter1 = 'MSG3-SEVI-MSGNDVE-0100-0100-201'
        filter2 = '.h5'
    print(downloadDir)
    folder = typeData
    print('--------------')
    #ftp.retrlines("LIST")
    ftp.cwd("eumetsat")
    ftp.cwd(folder)
    # ftp.cwd("subFolder")
    # or ftp.cwd("folderOne\\subFolder")
    print('--------------')
    #ftp.retrlines("LIST")
    listing = []
    lstFiles = []
    ftp.retrlines("LIST", listing.append)
    poslist = -1
    filename = []
    print('--------------')
    print('Number of files: ' + str(len(listing)))
    #print(listing)
    while(poslist < len(listing)):

        try:
            while True:
                poslist = poslist + 1
                print('.... File: ' + str(poslist))
                words = listing[poslist].split(None, 8)
                filename = words[-1].lstrip()
                # Filter to year data
                #print(filename)
                if(filter1 in str(filename) and  filter2 in str(filename)):
                    #print('!! pass in filters !!')
                    break
        except:
            if (poslist > len(listing)):
                printlog(True,'Finish copy files')
                break
        # download the file
        fileType = True
        #print(filename)
        try:
            id = filename.index('-201')
            #print('id in file')
        except:
            fileType = False
            printlog(True, 'File type different!')
        if(fileType == True):
            year  = filename[id+1:id+5]
            month = filename[id+5:id+7]
            subfolder = year + '\\' + month + '\\'
            #print(subfolder)
            if not(os.path.isdir(downloadDir + subfolder)):
                os.makedirs(downloadDir + subfolder)
            if(os.path.isfile(downloadDir + subfolder + filename) == True):
                printlog(True, 'File yet copy: ' + filename)
            else:
                printlog(True, 'File to copy: ' +  filename)
                lstFiles.append(filename)

    ftp.quit()
    return (lstFiles)

class myThread (threading.Thread):
    def __init__(self, threadID, name, q, typeData, downloadDir):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.q = q
        self.downloadDir = downloadDir
        self.typeData = typeData
    def run(self):
        printlog (True, "Starting " + self.name)
        #printlog(True, "q: " + str(self.q) + " typeData: " + self.typeData + " donwloadDir: " + self.downloadDir)
        process_data(self.name, self.q, self.typeData, self.downloadDir)
        printlog (True, "Exiting " + self.name)

def process_data(threadName, q, typeData, downloadDir):
    global queueLock
    global workQueue
    global exitFlag
    #printlog(True, "Enter in process_data. exitFlag: " + str(exitFlag))

    while not exitFlag:
        #printlog(True, "exitFlag is False")
        queueLock.acquire()
        #printlog(True, "\nqueueLock acquire")
        if not workQueue.empty():
            #printlog(True, "Queue is not empty")
            filename = q.get()
            queueLock.release()
            #printlog(True, "\nqueueLock release")
            printlog (True, "\n%s processing %s" % (threadName, filename))
            copyFTPfile(typeData, filename, downloadDir)
        else:
            #printlog(True, "Queue is empty")
            try:
                queueLock.release()
                #printlog(True, "\nqueueLock release")
            except:
                printlog(True, 'Fail to release queueLock !')
        time.sleep(1)
    printlog(True, "exitFlag is True")

def copyFTPfile(typeData, inputFile, outputDir):
    try:
        ftp = FTP("ftp.lapig.iesa.ufg.br", "lapig-ftp", "lapig123ftp")
    except:
        printlog(True, 'Error! Don''t have Internet Connection.')
        return(-1)
    try:
        ftp.login()
        printlog(False, 'Logged in ftp server')
        printlog(False, ftp.getwelcome())
    except:
        if(ftp.getwelcome()!= None):
            printlog(False, 'Logged in ftp server')
            printlog(False, ftp.getwelcome())
        else:
            printlog(False, 'Error! Verify user name and password')
            return(-2)
        printlog(False, 'Finally in login')
    folder = typeData
    ftp.cwd("eumetsat")
    ftp.cwd(folder)
    id = inputFile.index('-201')
    year  = inputFile[id+1:id+5]
    month = inputFile[id+5:id+7]
    subfolder = year + '\\' + month + '\\'
    outputFile = outputDir + subfolder + inputFile
    #printlog(True, "Create FTP data file")
    lf = open(outputFile, "wb")
    ftp.retrbinary("RETR " + inputFile, lf.write, 8*1024)
    lf.close()

class main:

    global exitFlag
    global threads
    global workQueue
    threadID = 1

    argDict = mapDict(argv, usage)

    if "-ty" in argDict and "-fo" in argDict:
        typeData = argDict["-ty"]
        downloadDir = argDict["-fo"]
    else:
        exit(usage)     
    if not(typeData == 'GRB' or typeData == 'HRIT' or typeData == 'HDF'):
        exit(usage)
    #print('-------- teste')
    #print(downloadDir)

    while(1):
        gc.collect()
        lst = ftpCopy(downloadDir, typeData)
        print('List files: ' + str(len(lst)))
        #print(lst)

        try:
            if(lst == -1):
                printlog(True, 'Waiting internet connection')
                sleepTime(0.1)
        except:
            break
        try:
            if(lst == []):    
                printlog(True, 'Waiting new files')
                sleepTime(1)
        except:
            break
        printlog(True, "----- Create new threads")
        for tName in threadList:
            thread = myThread(threadID, tName, workQueue, typeData, downloadDir)
            thread.start()
            threads.append(thread)
            threadID += 1

        printlog(True, "----- Fill the queue and wait case many files")
        
        x = 0
        while(x < len(lst)):
            queueLock.acquire()
            if not workQueue.full():
                filename =  lst[x]
                #printlog(True, 'x: ' + str(x) + '\tfilename: ' + filename)
                workQueue.put(filename)
                x = x + 1
                
            else:
                time.sleep(30)
            queueLock.release()

        printlog(True, "----- Wait for queue to empty")
        while not workQueue.empty():
            #print('Files in Queue: ' + str(workQueue.qsize()))
            #time.sleep(5)
            pass
        printlog(True, '----- Queue is empty')

        # Notify threads it's time to exit
        exitFlag = 1

    # Wait for all threads to complete
        for t in threads:
            t.join()
        printlog(True, "Exiting Main Thread")
        printlog(True, 'Waiting new files')
        sleepTime(1)
    
