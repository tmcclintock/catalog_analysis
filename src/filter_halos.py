"""
This is a stand-alone version of the filter_halos
function that takes in a halo file and produces
a new halo file without halos that aren't massive (Np >= 200)
enough that are not subhalos (pid > -1.0).
"""

def filter_halos(halo_path,filtered_halo_path):
    print "Filtering halo list."
    data = open(halo_path,"r")
    header = data.readline()
    outdata = open(filtered_halo_path,"w")
    outdata.write("#M X Y Z\n")
    for line in data:
        ID,DID,M,Vmax,Vrms,R200,Rs,Np,x,y,z,vx,vy,vz,pid = [float(item) for item in line.split()]
        if pid < 0.0 and Np >= 200:
            outdata.write("%e %e %e %e\n"%(M,x,y,z))
    outdata.close()
    print "\tHalos filtered."
    return
