export TrigChannel="0x22" #ch1=0x21, ch2=0x22,ch3=0x24,ch4=0x28,ext=0x30
export TrigLevel=100
export Offset=350 #delay before sampling begins
export Polarity="positive"
export SampleRate=2 #0=700MS/s, 1=1GS/s, 2=2GS/s, 3=5GS/s

rcdaq_client load librcdaqplugin_drs.so
rcdaq_client create_device device_drs -- 1 1001 $TrigChannel $TrigLevel $Polarity $Offset $SampleRate
rcdaq_control.pl &
