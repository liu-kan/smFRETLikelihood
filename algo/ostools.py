import platform,os

def hyperthreadingPerCore():
    sysstr = platform.system()
    thPerCpu=1
    if (sysstr == "Linux"):
        thPerCpu=int(os.popen("LC_ALL=C lscpu |grep Thread | awk '{print $4}'").readline().strip())
    elif (sysstr =="Windows"):
        from win32com.client import GetObject
        winmgmts_root = GetObject("winmgmts:root\cimv2")
        cpus = winmgmts_root.ExecQuery("Select * from Win32_Processor")
        for cpu in cpus:
            print('on "{}", hyperthreading is '.format(cpu.DeviceID), end='')
            if cpu.NumberOfCores < cpu.NumberOfLogicalProcessors:
                print('active')
                thPerCpu=cpu.NumberOfLogicalProcessors/cpu.NumberOfCores
            else:
                print('inactive')
    return thPerCpu
