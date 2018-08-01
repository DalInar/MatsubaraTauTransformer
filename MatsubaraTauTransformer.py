import numpy as np
import argparse

def G_tau_to_mat(taus, Gtau, beta, n_mat):
    mats = (2*np.arange(0, n_mat)+1)*np.pi/beta
    dt = taus[1]-taus[0]
    
    Gmat = np.zeros(n_mat, dtype=complex)  
    for i in range(0, n_mat):
        Gmat[i] = dt*np.sum(np.exp(1j*mats[i]*taus)*Gtau)
    
    return (mats, Gmat.real, Gmat.imag)

def G_mat_to_tau(mats, Gmat_real, Gmat_imag, beta, ntau):
    n_tau=ntau

    #Calculate the a1 coefficient for the high frequency tail, i.e. G ~ a1/(i\omega_n)
    a1 = -Gmat_imag[-1]*mats[-1]
    print("a1=",a1)

    taus = np.linspace(0.,beta,n_tau)
    n_mat = len(Gmat_real)
    iomega_n = 1j*(2*np.arange(0, n_mat)+1)*np.pi/beta
    Gtau = np.zeros(n_tau)
    Gmat = Gmat_real+1j*Gmat_imag
    for i in range(0, n_tau):
        tau = taus[i]
        Gtau[i] = (2./beta)*np.sum(np.exp(-iomega_n*tau)*(Gmat - a1/iomega_n)).real
      
    return (taus, Gtau-a1/2.)

def main():
	parser = argparse.ArgumentParser(description="Hi! I transform greens functions between matsubara and tau space.")
	parser.add_argument('input_file', help='File to read in data from')
	parser.add_argument('output_file', help='File to output results to')
	parser.add_argument('input_fnc_type', help="'greens' or 'selfenergy'")
	parser.add_argument('input_dataspace', help="'matsubara' (freq, ReG, ImG) or 'tau' (tau, ReG)")
	parser.add_argument('beta', help="beta = inverse temperature")
	args = parser.parse_args()

	beta = float(args.beta)
	if(args.input_dataspace == 'matsubara'):
		print("Performing transform from matsubara to tau space.")
		G_omega = np.loadtxt(args.input_file)
		mats = G_omega[:,0]
		G_omega_real = G_omega[:,1]
		G_omega_imag = G_omega[:,2]
		beta_calc = np.pi/mats[0]
		
		if(abs(beta-beta_calc)>0.01):
			print('Command line beta ',beta,' does not match the beta calculated from the matsubara frequencies ',beta_calc)
			exit()
		constant_shift = 0
		if(args.input_fnc_type == 'selfenergy'):
			constant_shift = G_omega_real[-1]
			G_omega_real -= constant_shift
			print("Constant term in Re(selfenergy) = ", constant_shift)
		
		(taus, G_tau_real) = G_mat_to_tau(mats, G_omega_real, G_omega_imag, beta, len(mats)*2)
		G_tau = np.zeros(shape=(len(taus), 2))
		G_tau[:,0] = taus
		G_tau[:,1] = G_tau_real
		np.savetxt(args.output_file, G_tau)
		print("Done!")

	elif(args.input_dataspace == 'tau'):
		print("Performing transform from tau to matsubara.")
		G_tau = np.loadtxt(args.input_file)
		n_mat = 512
		(mats, G_omega_real, G_omega_imag) = G_tau_to_mat(G_tau[:,0], G_tau[:,1], beta, n_mat)
		G_omega = np.zeros(shape=(len(mats),3))
		G_omega[:,0] = mats
		G_omega[:,1] = G_omega_real
		G_omega[:,2] = G_omega_imag
		np.savetxt(args.output_file, G_omega)
		print("Done!")

	else:
		print("Invalid 'input_dataspace', must be either 'matsubara' or 'tau'.")


if __name__=="__main__":
	main()

