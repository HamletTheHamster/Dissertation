import csv
import visa
import time
import datetime
import numpy as np
import matplotlib.pyplot as plt
import zhinst.ziPython, zhinst.utils
import os
import winsound

daq = zhinst.ziPython.ziDAQServer('localhost', 8005)
device = zhinst.utils.autoDetect(daq)

startTime = datetime.datetime.now()

dwell = .1
F_lockin = 45.0000E6
F_AOM = 40.000800e6
chan = 1
ch   = str(chan-1)
rate=2e9 # rate, per millisecond
RBW = 100 # lock-in bandwidth
tc=1/(2*np.pi*RBW) # time constant
lockRange = .04
order = 8

exp_set = [
       [['/', device, '/sigins/',ch,'/diff'], 0],
       [['/', device, '/sigins/',ch,'/imp50'], 1],
       [['/', device, '/sigins/',ch,'/ac'], 1],

       # Want range as low as possible without clipping data
       # (red Over light on Lock-in box)
       [['/', device, '/sigins/',ch,'/range'], lockRange],
       [['/', device, '/demods/',ch,'/order'], order],
       [['/', device, '/demods/',ch,'/timeconstant'], tc/3.33],
       [['/', device, '/demods/',ch,'/rate'], rate],
       [['/', device, '/demods/',ch,'/adcselect'], chan-1],
       [['/', device, '/demods/',ch,'/oscselect'], chan-1],
       [['/', device, '/demods/',ch,'/harmonic'], 1],
       [['/', device, '/oscs/',ch,'/freq'], F_lockin],
   ]

daq.set(exp_set);
time.sleep(.001)

path = '/%s/demods/%d/sample' % (device, 0) # (device, demod_index)
daq.subscribe(path)
daq.flush()
daq.sync()

def measureSynchronousFeedback(daq, device, channel, frequency):
   c=str(channel-1) #return a string of an object. =channel-1=-1 for channel 0

# Poll the subscribed data from the data server. Poll will block and record
# for poll_length seconds.
    daq.sync()
    daq.flush()
    poll_length = dwell  # [s]
    poll_timeout = 10  # [ms]
    poll_flags = 0
    poll_return_flat_dict = True
    data = daq.poll(
        poll_length,
        poll_timeout,
        poll_flags,
        poll_return_flat_dict
    )
    assert data, """poll() returned an empty data dictionary,
                    did you subscribe to any paths?"""

    # Access the demodulator sample using the node's path.
    sample = data[path] # Defines the sample as the data from defined path

    # Calculate the demodulator's magnitude and adds it to the dict.
    global sampleR
    sampleR = np.abs(sample['x'] + 1j*sample['y']) #y-axis #magnitude

    clockbase = float(daq.getInt('/%s/clockbase' % device))

    # Convert timestamps from ticks to seconds via clockbase.
    t = (sample['timestamp'] - sample['timestamp'][0])/clockbase

resources = visa.ResourceManager()
resources.list_resources()

#open the signal generator
#gen=resources.open_resource('USB0::0x03EB::0xAFFF::6C2-0A2A2000A-0374::INSTR')
gen=resources.open_resource('USB0::0x03EB::0xAFFF::6C2-0A2B2000A-0430::INSTR')

gen.write('FREQUENCY:MODE CW DUAL')
time.sleep(0.01)
gen.write('OUTPUT1 1')
time.sleep(0.01)
gen.write('OUTPUT2 1')

# run a cycle of frequencies for each output
Fstart = 9.0*1E9
Fstop  = 9.28*1E9
F_step = 0.00500456*1E9

N_steps = np.int((Fstop-Fstart)/F_step)

NoRealizations = 5
dataR = [0]*N_steps # demod data
dataF = [0]*N_steps # Freq. axis
stdDevOfMeanR = [0]*N_steps # to hold standard deviation of dwell-time data
f2 = Fstart
F = Fstart


run = 1
folderName = startTime.strftime("%y-%m-%d ") + """1cm UHNA3"""
while os.path.exists(folderName + "/" + str(run) + "/signal.csv"):
    run += 1
folderRunName = folderName + "/" + str(run)

signalData = 1
takeBkrdData = 0

aveDataR = [0]*N_steps

if not os.path.exists(folderRunName):
    os.makedirs(folderRunName)

if not os.path.exists(folderName + "/meta.csv"):
    meta = open(folderName + "/meta.csv", "a")
    meta.write("""Date,Label,Sets,Start Time,End Time (hr:min:sec),Pump,Stokes,
                Probe,Frequency,Signal,Range,Dwell,Bandwidth,Data Rate,Order,
                Start Frequency,Stop Frequency,Step,Num Avgs,Pump Laser,Probe
                Laser,Probe Filter,Stokes Filter,Notes\n""")
    meta.close()

if takeBkrdData == 1:
    bkrdStartTime = datetime.datetime.now()

    bgDataR = np.zeros((NoRealizations, N_steps))

for kk in range(0,NoRealizations):
    F = Fstart
    print(kk)
    for jj in range(0,N_steps):
        f1 = F
        f2 = F + F_AOM - F_lockin
        time.sleep(0.1)
        gen.write("SOURce1: FREQuency:CW "+ str(f1)) # Hz
        time.sleep(1e-3)
        gen.write("SOURce2:FREQuency:CW "+ str(f2)) # Hz
        time.sleep(1e-3)
        measureSynchronousFeedback(daq, device, 1, F_lockin)
        time.sleep(1e-3)

        F=F+F_step

        dataF[jj] = Fstart+jj*F_step
        ff1 = f1*10**-6
        ff2 = f2*10**-6

        dataR[jj]=np.mean(sampleR)
        stdDevOfMeanR[jj] = np.std(sampleR)/np.sqrt(len(sampleR))

    #-----------------save run data (dataR and stdDevR)------------------------#
    if signalData == 1:
        runSigDir = folderRunName + "/Runs/Signal/"
        if not os.path.exists(runSigDir):
            os.makedirs(runSigDir)

        csvfile = runSigDir + "Run " + str(kk) + ".csv"
        with open(csvfile, "w") as output:
            writer = csv.writer(output, delimiter=',')
            writer.writerow(['Sig','Std Dev'])
            for sig, std in zip(dataR, stdDevOfMeanR):
                writer.writerow([sig, std])

    elif takeBkrdData == 1:
        runBgDir = folderRunName + "/Runs/Background/"
        if not os.path.exists(runBgDir):
            os.makedirs(runBgDir)

        csvfile = runBgDir + "Run " + str(kk) + ".csv"
        with open(csvfile, "w") as output:
            writer = csv.writer(output, delimiter=',')
            writer.writerow(['Sig','Std Dev'])
            for sig, std in zip(dataR, stdDevOfMeanR):
                writer.writerow([sig, std])
    #--------------------------------------------------------------------------#

    aveDataR = (np.array(dataR)+kk*np.array(aveDataR))/(kk+1)

    if takeBkrdData == 1:

        bgDataR[kk] = dataR
        background = aveDataR

        plt.figure()
        plt.grid(True)
        plt.plot(dataF, background)
        plt.title('Background Demodulator data')
        plt.show()

    if takeBkrdData == 0:
        plt.figure()
        plt.grid(True)
        plt.plot(dataF, aveDataR-background)
        plt.title('Signal - Background Demodulator data')
        plt.xlabel('F')
        plt.ylabel('R')
        plt.show()

    #-----------------save subtracted run data (dataR)-------------------------#
    if signalData == 1:
        runSubtrDir = folderRunName + "/Runs/Subtracted/"
        if not os.path.exists(runSubtrDir):
            os.makedirs(runSubtrDir)

        csvfile = runSubtrDir + "Run " + str(kk) + ".csv"
        with open(csvfile, "w") as output:
            writer = csv.writer(output, lineterminator='\n')
            for val in dataR - bgDataR[kk]:
                writer.writerow([val])
    #--------------------------------------------------------------------------#

if signalData == 1:
    csvfile=folderRunName+"/signal.csv"
    with open(csvfile, "w") as output:
        writer=csv.writer(output, lineterminator='\n')
        for val in aveDataR-background:
            writer.writerow([val])

    csvfile=folderRunName+"/frequency.csv"
    with open(csvfile, "w") as output:
        writer=csv.writer(output, lineterminator='\n')
        for val in dataF:
            writer.writerow([val])

if not os.path.exists(folderRunName + "/timestamp.csv"):
    timestampf = open(folderRunName + "/timestamp.csv", "a")
    timestampf.write("""Date,Label,Run,Data,Start Time,End Time (hr:min:sec),
                    Pump,Stokes,Probe,Frequency,Signal,Range,Dwell,Bandwidth,
                    Data Rate,Order,Start Frequency, Stop Frequency,Step,Num
                    Avgs,Notes\n""")
    timestampf.close()

if takeBkrdData == 1:
    bkrdEndTime = datetime.datetime.now()
    timestampf = open(folderRunName + "/timestamp.csv", "a")
    timestampf.write(
        bkrdStartTime.strftime("%y-%m-%d,,") + str(run) + ",Background" +
        bkrdStartTime.strftime(",%H:%M:%S") +
        bkrdEndTime.strftime(",%H:%M:%S,,,,,,")
    )
    timestampf.write(
        str(lockRange) + "," + str(dwell) + "," + str(RBW) + "," +
        str("{:e}".format(rate)) + "," + str(order) + "," +
        str("{:e}".format(Fstart)) + "," + str("{:e}".format(Fstop)) + "," +
        str("{:e}".format(F_step)) + "," + str(NoRealizations) +",\n"
    )
    timestampf.close()

if signalData == 1:
    endTime = datetime.datetime.now()
    timestampf = open(folderRunName + "/timestamp.csv", "a")
    timestampf.write(
        startTime.strftime("%y-%m-%d,,") + str(run) + ",Signal" +
        startTime.strftime(",%H:%M:%S") + endTime.strftime(",%H:%M:%S,,,,,")
    )
    timestampf.write(str(run) + "/frequency.csv," + str(run) + "/signal.csv,")
    timestampf.write(
        str(lockRange) + "," + str(dwell) + "," + str(RBW) + "," +
        str("{:e}".format(rate)) + "," + str(order) + "," +
        str("{:e}".format(Fstart)) + "," + str("{:e}".format(Fstop)) + "," +
        str("{:e}".format(F_step)) +  "," + str(NoRealizations) +",\n"
    )
    timestampf.close()

    meta = open(folderName + "/meta.csv", "a")
    meta.write(
        startTime.strftime("%y-%m-%d,,") + str(run) +
        startTime.strftime(",%H:%M:%S") + endTime.strftime(",%H:%M:%S,,,,")
    )
    meta.write(str(run) + "/frequency.csv," + str(run) + "/signal.csv,")
    meta.write(
        str(lockRange) + "," + str(dwell) + "," + str(RBW) + "," +
        str("{:e}".format(rate)) + "," + str(order) + "," +
        str("{:e}".format(Fstart)) + "," + str("{:e}".format(Fstop)) + "," +
        str("{:e}".format(F_step)) +  "," + str(NoRealizations) +",,,,,\n"
    )
    meta.close()
