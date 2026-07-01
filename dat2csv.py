import time, sys, os

def print_buffer(buff):
   pt = 0
   for i in range(0,416):
      if pt == 0:
         print("%03X::%02X " %(i,buff[i]), end="")
      else:
         print("%02X " %(buff[i]), end="")
      pt = pt+1
      if pt == 16:
         print("")
         pt = 0

def extract_data_2_csv(datfilename, csvfilename):
   err_frameNb = 0
   filesize = os.stat(datfilename).st_size
   NbTrig = 0;
   oldtimestamp = -1;
   ffirst = True;
   firstevent = 0
   lastevent = 0
   min_ecart = -1;
   max_ecart = 0;
   if filesize >= 416:
      if csvfilename.upper() != "NONE":
         outFile = open(csvfilename,'w')
         txt = "Conv  {:>16}".format(os.path.basename(datfilename))
         #print("Convert {} to {}".format(os.path.basename(datfilename),csvfilename), end = "")
      else:
         txt = "Check {:>16}".format(os.path.basename(datfilename))
         #print("Check {:>20}".format(os.path.basename(datfilename)))
      inFile = open(datfilename, 'rb')
      buff = inFile.read(416)
      packet = 0
      while(buff):
         packet = packet + 1
         #print("Packet N°:%d" %(packet))
         metaType = int.from_bytes(buff[0:4],byteorder='little',signed=False)
         dataSource = int.from_bytes(buff[4:8],byteorder='little',signed=False)
         frameNb = (dataSource >> 8) & 0xFFFF
         dataSource = dataSource & 0xFF0000FF
         #print('metaType = 0x%08X, dataSource = 0x%08X' %(metaType,dataSource))
         if (metaType == 0x000068C2) and (dataSource == 0x01000001):
            NbTrig = NbTrig + 1
            #print("Nb Trame Decoded = %15d\r"%(NbTrig), end="")
            ts_lsb = int.from_bytes(buff[8:12],byteorder='little',signed=False)
            ts_msb = int.from_bytes(buff[12:16],byteorder='little',signed=False)
            #timestamp = get_timestamp(ts_lsb,ts_msb)
            timestamp = get_timestamp(ts_lsb,ts_msb)
            if ffirst == False:
               if frameNb != ((lastevent + 1) & 0xFFFF):
                  err_frameNb = err_frameNb + 1;
               lastevent = frameNb
               if timestamp >= oldtimestamp:
                  ecart = timestamp - oldtimestamp;
               else:
                  print("\n{}-{} -> {} : {}".format(ts_msb,ts_lsb,timestamp,oldtimestamp))
                  ecart = (0x10000000000000000 - oldtimestamp)+timestamp
               if ecart > max_ecart:
                  max_ecart = ecart
               if ecart < min_ecart or min_ecart == -1:
                  min_ecart = ecart
            else:
               firstevent = frameNb
               lastevent = frameNb
            
            ffirst = False
            oldtimestamp = timestamp
            for i in range(0,16):
               header = buff[(i*25)+16]
               if header != 0xAC:
                  print(txt)
                  print("ERROR::Asic frame header (index = 0x%03X : %02X %02X)" %(((i*25)+16),header,buff[(i*25)+17]))
                  print_buffer(buff)
                  return 0
               elif csvfilename.upper() != "NONE":
                  #print("Asic frame header (index = %d : 0x%02X)" %(((i*25)+16),header))
                  txtcsv = str(timestamp) + ";"
                  for j in range(0,8):
                     index = (i*25)+17+(j*3)
                     data = int.from_bytes(buff[index:index+3],byteorder='big',signed=False)
                     #data = int.from_bytes(buffer[index:index+3],byteorder='little',signed=False)
                     # Extract Bits from data
                     str_data = format(data, '024b')[::-1]
                     # New format
                     tdc = int(str_data[0:12][::-1],2)  # Bit-11:0
                     hit = int(str_data[12],2)          # Bit-12
                     adc = int(str_data[13:21][::-1],2) # Bit-20:13
                     if j != 7:
                        txtcsv = txtcsv + ("%d;%d;%d;" %(tdc,adc,hit))
                     else:
                        txtcsv = txtcsv + ("%d;%d;%d" %(tdc,adc,hit))
                  outFile.writelines(txtcsv)
                  outFile.write('\n')            
            #print("FRAME Nb = {}".format(frameNb))
         else:
            print(txt)
            print("%d -> BAD HEADER 0x%08X 0x%08X"%(frameNb,metaType,dataSource))
         oldbuff = buff
         buff = inFile.read(416)
      if err_frameNb > 0:
         print("{} | {:0>8} | {:>5} -> {:>5}*| {} :: {} ".format(txt,packet,firstevent,lastevent,min_ecart,max_ecart))
      else:
         print("{} | {:0>8} | {:>5} -> {:>5} | {} :: {} ".format(txt,packet,firstevent,lastevent,min_ecart,max_ecart))
      return packet
# =======================================================
#
# =======================================================
def get_timestamp(lsb, msb):
   num = (msb << 32) + lsb
   return num
# =======================================================
#
# =======================================================
def get_timestamp_gray(lsb, msb):
   num = (msb << 32) + lsb
   ures = num ^ (num >> 32)
   ures ^= (ures >> 16)
   ures ^= (ures >> 8)
   ures ^= (ures >> 4)
   ures ^= (ures >> 2)
   ures ^= (ures >> 1)
   return ures
# =======================================================
#
# =======================================================
if __name__ == "__main__":
   total_packet = 0
   print("-----------------------------------------------------------------------")
   print("      File Names       | Nb Trigs |  Event Number  | TimeStamp Min-Max")
   print("-----------------------|----------|----------------|-------------------")
   if len(sys.argv) == 3:
      datfile = sys.argv[1]
      csvfile = sys.argv[2]
      ret = extract_data_2_csv(datfile, csvfile)
      if ret > 0:
         total_packet = total_packet + ret
   elif len(sys.argv) > 4:
      datfile = sys.argv[1]
      csvfile = sys.argv[2]
      begin_index = int(sys.argv[3])
      end_index = int(sys.argv[4])
      for i in range(begin_index,end_index+1):
         datf = datfile + ("_%d.dat" %(i))
         if csvfile.upper() != "NONE":
            csvf = csvfile + ("_%d.csv" %(i))
         else:
            csvf = "NONE"
         nbTry = 0
         while nbTry < 30:
            if os.path.isfile(datf):
               ret = extract_data_2_csv(datf, csvf)
               if ret > 0:
                  total_packet = total_packet + ret
               nbTry = 65535
            else:
               nbTry = nbTry + 1
               time.sleep(1)
         if nbTry != 65535:
            print("File {} not exist".format(datf))
            break

   print("TOTAL TRIGGER(s) = {}".format(total_packet))
