package hetmm;

import java.util.List;
import java.util.Random;

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.inference.State;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;


@Description("Likelihood for fitting a cocktail of enzyme kinetics")
public class HeterogeneousLikelihood extends HomogeneousLikelihood {

	final public Input<RealParameter> vmaxVarInput = new Input<>("vmaxVar", "variance of vmax", Input.Validate.REQUIRED);
	final public Input<RealParameter> kmVarInput = new Input<>("kmVar", "variance of Km", Input.Validate.REQUIRED);
	final public Input<RealParameter> rhoInput = new Input<>("rho", "correlation between Km and vmax", Input.Validate.REQUIRED);
	final public Input<Integer> nStepsInput = new Input<>("nsteps", "grid size in numerical integration", 100);
	
	
	final public Input<BooleanParameter> indicatorInput = new Input<>("indicator", "indicator is true if this model is being used, false for HomogeneousLikelihood");
	
	
	
	final double NUMBER_OF_SD_IN_GRID = 3;
	
	
	MultivariateNormalDistribution mvn;
	double[] vmaxGridRange = new double[2];
	double[] kmridRange = new double[2];
	double[][] gridDensity;
	
	int ngrid = 0;
	
	
	@Override
	public void initAndValidate() {
		
		super.initAndValidate();
		
		this.ngrid = nStepsInput.get();
		gridDensity = new double[ngrid+1][ngrid+1];
		
		rhoInput.get().setBounds(-1.0, 1.0);
		
		prepare();
		
		
		
		
		//throw new IllegalArgumentException("Please ensure a is >= 0");
		
		

	}
	
	
	@Override
	public boolean prepare() {
		
		if (indicatorInput.get() != null && !indicatorInput.get().getValue()) {
			return super.prepare();
		}
		
		// Parameters
		double kmMean = kmInput.get().getValue();
		double vmaxMean = vmaxInput.get().getValue();
		double kmVar = kmVarInput.get().getValue();
		double vmaxVar = vmaxVarInput.get().getValue();
		double rho = rhoInput.get().getValue();
        
        double covar =	rho* Math.sqrt(kmVar)*Math.sqrt(vmaxVar);
		

		
		// Transform means into log space
		kmMean = Math.log(kmMean) - kmVar/2;
		vmaxMean = Math.log(vmaxMean) - vmaxVar/2;
		

		
		if (Math.abs(covar) > kmVar*vmaxVar) {
			Log.warning("illegal covariance " + kmMean + " " + vmaxMean + " " + kmVar + " " + vmaxVar + " " + covar);
			return false;
		}
		
		//
		
		
		// Mean vector
		double[] means = new double[] { kmMean, vmaxMean };
		
		// Covariance matrix
		double[][] variances = new double[2][2];
		variances[0][0] = kmVar;
		variances[0][1] = covar;
		variances[1][0] = covar;
		variances[1][1] = vmaxVar;
		
		//Log.warning(" covariance " + kmMean + " " + vmaxMean + " " + kmVar + " " + vmaxVar + " " + covar);
		this.mvn = new MultivariateNormalDistribution(means, variances);
		
		
		
		// Grid
		vmaxGridRange[0] = vmaxMean - Math.sqrt(vmaxVar)*NUMBER_OF_SD_IN_GRID;
		vmaxGridRange[1] = vmaxMean + Math.sqrt(vmaxVar)*NUMBER_OF_SD_IN_GRID;
		
		kmridRange[0] = kmMean - Math.sqrt(kmVar)*NUMBER_OF_SD_IN_GRID;
		kmridRange[1] = kmMean + Math.sqrt(kmVar)*NUMBER_OF_SD_IN_GRID;
		

		// Density
		double vmaxdt = (vmaxGridRange[1]-vmaxGridRange[0])/this.ngrid;
		double kmdt = (kmridRange[1]-kmridRange[0])/this.ngrid;;
		for (int i = 1; i <= this.ngrid; i ++) {
			
			double vmaxLower = vmaxGridRange[0] + (i-1) * vmaxdt;
			double vmaxUpper = vmaxGridRange[0] + (i) * vmaxdt;
			double vmax = (vmaxUpper+vmaxLower)/2;
			
			// Transform back into positive domain
			vmax = Math.exp(vmax);
			
			for (int j = 1; j <= this.ngrid; j ++) {
				
				double kmLower = kmridRange[0] + (j-1) * kmdt;
				double kmUpper = kmridRange[0] + (j) * kmdt;
				double km = (kmUpper+kmLower)/2;
				
				
				// Transform back into positive domain
				km = Math.exp(km);
				
				// Density
				double p1 = this.mvn.density(new double[] { kmLower, vmaxLower } );
				double p2 = this.mvn.density(new double[] { kmUpper, vmaxUpper } );
				double p = (vmaxUpper-vmaxLower) * (kmUpper-kmLower) * (p2+p1)/2;
				gridDensity[i][j] = p;
				
			}
			
		}
		
		
		return true;
		
	}
	
	/**
	 * The expected rate velocity under the model, as a function of reactant concentration
	 * @param a
	 * @return
	 */
	public double getExpectedV(double a) {
		
		
		if (indicatorInput.get() != null && !indicatorInput.get().getValue()) {
			return super.getExpectedV(a);
		}
		
		
		double v = 0;
		
		double hill = getHill();
		if (hill != 1) {
			a = Math.pow(a, hill);
		}
		
		double pSum = 0;
		
		double vmaxdt = (vmaxGridRange[1]-vmaxGridRange[0])/this.ngrid;
		double kmdt = (kmridRange[1]-kmridRange[0])/this.ngrid;;
		
		// Find grid boundaries
		for (int i = 1; i <= this.ngrid; i ++) {
			
			double vmaxLower = vmaxGridRange[0] + (i-1) * vmaxdt;
			double vmaxUpper = vmaxGridRange[0] + (i) * vmaxdt;
			double vmax = (vmaxUpper+vmaxLower)/2;
			
			
			
			// Transform back into positive domain
			vmax = Math.exp(vmax);
			
			for (int j = 1; j <= this.ngrid; j ++) {
				
				
				
				
				double kmLower = kmridRange[0] + (j-1) * kmdt;
				double kmUpper = kmridRange[0] + (j) * kmdt;
				double km = (kmUpper+kmLower)/2;
	
				
				// Transform back into positive domain
				km = Math.exp(km);
				
				if (hill != 1) {
					km = Math.pow(km, hill);
				}
				
				
				// Probability density
				double p = gridDensity[i][j];
				
				
				pSum += p;
				
				// Calculate v
				double vThis = (vmax*a) / (km + a);
				
				//Log.warning("\tp= " + p + " v(a) = v(" + a + ")=" + vThis + " ------ >" + (vThis*p));
				
				v += vThis*p;
				
				
			}
			
			
		}
		
		
		
		v = v / pSum;
		
		//Log.warning("psum = " + pSum + " v(a) = v(" + a + ")=" + v + " | " + kmridRange[0]);
		
		//this.mvn.density(null)
		
		// Integrate across multivariate normal distribution
		//this.mvn.
		
		//return (vmax*a) / (km + a);
		
		return v;
	}
	
	
	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub
		
	}

}
