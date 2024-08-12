import numpy as np
import datetime as dt
import h5py
from apexpy import Apex
import matplotlib.pyplot as plt
from loess.loess_1d import loess_1d

def main():
    def get_h5_data(var):
        data = np.array(h5f.get(var)).flatten()
        data = data[track_id-window:track_id+window]
        return data

    path = 'data/swop/swarm_data/'
    sims = open('data/swop/sim_data.txt','r').readlines()

    for sim in sims:
        sim_data = sim.split()
        filename = sim_data[0]
        print(filename)
        sim_time = np.datetime64(sim_data[1])
        vx_choice = sim_data[2]

        h5f = h5py.File(path + filename + '.h5','r')
        track_time = np.array(h5f.get('Timestamp')).astype('datetime64[s]').flatten()
        track_id = np.argmin(np.abs(track_time-sim_time))
        window = 2*6*60 # 2 Hz x 6 minutes

        track_glat = get_h5_data('GeodeticLatitude')
        track_glon = get_h5_data('Longitude') + 360
        track_alt = get_h5_data('GeodeticAltitude') / 1e3
        track_vxh = get_h5_data('Vixh')
        track_vxv = get_h5_data('Vixv')
        track_vy = get_h5_data('Viy')
        track_vz = get_h5_data('Viz')
        track_nE = get_h5_data('VsatE')
        track_nN = get_h5_data('VsatN')
        # track_nU = -get_h5_data('VsatC') # assuming nU = 0, ~0.1% of nN
        h5f.close()

        track_nnorm = np.sqrt(track_nE**2 + track_nN**2)
        track_nE = track_nE / track_nnorm
        track_nN = track_nN / track_nnorm
        # track_nU = track_nU / track_nnorm

        A = Apex(datetime64_to_datetime(sim_time))
        # track_alt_new = np.ones(track_alt.size)*110
        track_mlat, track_mlon = A.convert(track_glat,track_glon,'geo','apex',height=track_alt)
        track_110 = A.map_to_height(track_glat,track_glon,track_alt,110)
        track_glat_110 = track_110[0]
        track_glon_110 = track_110[1]
        track_mlon += 360

        loess_n = 15
        if vx_choice == 'h':
            _, track_vx_s, _ = loess_1d(track_mlat, track_vxh, npoints=loess_n)
        elif vx_choice == 'v':
            _, track_vx_s, _  = loess_1d(track_mlat, track_vxv, npoints=loess_n)
        elif vx_choice == 'b':
            track_vx = np.hstack((track_vxh,track_vxv))
            dummy_vx = np.hstack((track_mlat,track_mlat))
            _, track_vx_s, _  = loess_1d(dummy_vx, track_vx, npoints=2*loess_n, xnew=track_mlat)
        _, track_vy_s, _  = loess_1d(track_mlat, track_vy, npoints=loess_n)
        _, track_vz_s, _  = loess_1d(track_mlat, track_vz, npoints=loess_n)

        xlim0 = track_mlat[int(len(track_mlat)/2-120)]
        xlim1 = track_mlat[int(len(track_mlat)/2+120)]
        xlims = [xlim0, xlim1]

        if (vx_choice == 'h') | (vx_choice == 'b'):
            plt.plot(track_mlat,track_vxh, ':r', linewidth=0.5, label='Vixh')
        if (vx_choice == 'v') | (vx_choice == 'b'):
            plt.plot(track_mlat,track_vxv, ':m', linewidth=0.5, label='Vixv')
        plt.plot(track_mlat,track_vx_s, 'b', linewidth=0.5, label='Vix loess')
        plt.xlim(xlims)
        plt.ylim(-3e3,3e3)
        plt.xlabel('Magnetic latitude (°)')
        plt.ylabel('Ion flow (m/s)')
        plt.title(filename)
        plt.legend()
        plt.savefig(path + 'plots/' + filename + '_vx.png', dpi=600)
        plt.close()

        plt.plot(track_mlat,track_vy, ':m', linewidth=0.5, label='Viy')
        plt.plot(track_mlat,track_vy_s, 'r', linewidth=0.5, label='Viy loess')
        plt.plot(track_mlat,track_vz, ':c', linewidth=0.5, label='Viz')
        plt.plot(track_mlat,track_vz_s, 'b', linewidth=0.5, label='Viz loess')
        plt.xlim(xlims)
        plt.ylim(-3e3,3e3)
        plt.xlabel('Magnetic latitude (°)')
        plt.ylabel('Ion flow (m/s)')
        plt.title(filename)
        plt.legend()
        plt.savefig(path + 'plots/' + filename + '_vyvz.png', dpi=600)
        plt.close()

        track_vE = track_vx_s*track_nE + track_vy_s*track_nN
        track_vN = track_vx_s*track_nN - track_vy_s*track_nE
        track_vU = -track_vz_s
        track_vENU = np.vstack((track_vE,track_vN,track_vU))

        *_, e1, _, e3 = A.basevectors_apex(track_glat, track_glon, track_alt)
        x1 = -e3 / np.linalg.norm(e3, axis=0, keepdims=True)
        u2 = e1 - dot(x1,e1)*x1
        x2 = u2 / np.linalg.norm(u2, axis=0, keepdims=True)
        x3 = np.cross(x1.T, x2.T).T

        track_vmE = dot(track_vENU,x2)
        track_vmN = dot(track_vENU,x3)
        track_vmU = dot(track_vENU,x1)

        plt.plot(track_mlat,track_vE, ':m', linewidth=0.5, label='v_E geog')
        plt.plot(track_mlat,track_vmE, 'r', linewidth=0.5, label='v_E geom')
        plt.xlim(xlims)
        plt.ylim(-3e3,3e3)
        plt.xlabel('Magnetic latitude (°)')
        plt.ylabel('Ion flow (m/s)')
        plt.title(filename)
        plt.legend()
        plt.savefig(path + 'plots/' + filename + '_magE.png', dpi=600)
        plt.close()
        
        plt.plot(track_mlat,track_vN, ':m', linewidth=0.5, label='v_N geog')
        plt.plot(track_mlat,track_vmN, 'r', linewidth=0.5, label='v_N geom')
        plt.xlim(xlims)
        plt.ylim(-3e3,3e3)
        plt.xlabel('Magnetic latitude (°)')
        plt.ylabel('Ion flow (m/s)')
        plt.title(filename)
        plt.legend()
        plt.savefig(path + 'plots/' + filename + '_magN.png', dpi=600)
        plt.close()
        
        plt.plot(track_mlat,track_vU, ':m', linewidth=0.5, label='v_U geog')
        plt.plot(track_mlat,track_vmU, 'r', linewidth=0.5, label='v_U geom')
        plt.xlim(xlims)
        plt.ylim(-3e3,3e3)
        plt.xlabel('Magnetic latitude (°)')
        plt.ylabel('Ion flow (m/s)')
        plt.title(filename)
        plt.legend()
        plt.savefig(path + 'plots/' + filename + '.png', dpi=600)
        plt.close()

        h5f = h5py.File(path + filename + '.h5','a')
        write_h5_data(h5f,'MagneticLongitude',track_mlon)
        write_h5_data(h5f,'MagneticLatitude',track_mlat)
        write_h5_data(h5f,'ViE',track_vE)
        write_h5_data(h5f,'ViN',track_vN)
        write_h5_data(h5f,'ViU',track_vU)
        write_h5_data(h5f,'ViMagE',track_vmE)
        write_h5_data(h5f,'ViMagN',track_vmN)
        write_h5_data(h5f,'ViMagU',track_vmU)
        write_h5_data(h5f,'ids',np.array((track_id-window,track_id+window)))
        write_h5_data(h5f,'Longitude110km',track_glon_110)
        write_h5_data(h5f,'GeodeticLatitude110km',track_glat_110)
        h5f.close()

def dot(a1,a2):
    return np.sum(a1*a2,axis=0)

def datetime64_to_datetime(date):
    timestamp = (date - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')
    return dt.datetime.fromtimestamp(timestamp,tz=dt.timezone.utc)

def write_h5_data(file,name,data):
    try:
        file.create_dataset(name, data=data)
    except:
        print('{} already exists in {}.'.format(name,file))

if __name__ == '__main__':
    main()