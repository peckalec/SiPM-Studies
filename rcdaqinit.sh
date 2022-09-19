export TrigChannel="0x320" #ch1=0x21, ch2=0x22,ch3=0x24,ch4=0x28,ext=0x30,ch1&ch2=0x320
export TrigLevel=15
export Offset=100 #delay before sampling begins, 350 for 2 GS/s
export Polarity="positive"
export SampleRate=3 #0=700MS/s, 1=1GS/s, 2=2GS/s, 3=5GS/s
export Baseline=-450 #baseline shift


rcdaq_client load librcdaqplugin_drs.so
rcdaq_client create_device device_drs -- 1 1001 $TrigChannel $TrigLevel $Polarity $Offset $SampleRate 0 1024 $Baseline
rcdaq_control.pl &
