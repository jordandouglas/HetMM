package hetmm;

import org.apache.commons.math3.distribution.NormalDistribution;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;

@Description("Likelihood for fitting a mixture of enzymes with varying Km")
public class HKLikelihood extends HomogeneousLikelihood {
	
	
	final public Input<RealParameter> kmSDInput = new Input<>("kmSD", "SD of Km", Input.Validate.REQUIRED);
	final public Input<Integer> nStepsInput = new Input<>("nsteps", "grid size in numerical integration", 1000);
	final public Input<BooleanParameter> indicatorInput = new Input<>("indicator", "indicator is true if this model is being used, false for HomogeneousLikelihood");
	
	
	final double NUMBER_OF_SD_IN_GRID = 8;
	NormalDistribution dist;
	double[] kmridRange = new double[2];
	double[] gridDensity;
	
	int ngrid = 0;
	
	
	@Override
	public void initAndValidate() {
		
		super.initAndValidate();
		
		this.ngrid = nStepsInput.get();
		gridDensity = new double[ngrid+1];
		
		
		prepare();
		

	}
	
	

	@Override
	public boolean prepare() {
		
		if (indicatorInput.get() != null && !indicatorInput.get().getValue()) {
			return super.prepare();
		}
		
		// Parameters
		double kmMean = kmInput.get().getValue();
		double kmSD = kmSDInput.get().getValue();
        
		
		// Transform mean into log space
		kmMean = Math.log(kmMean) - kmSD/2;
		

		
		
		//Log.warning(" covariance " + kmMean + " " + vmaxMean + " " + kmVar + " " + vmaxVar + " " + covar);
		this.dist = new NormalDistribution(kmMean, kmSD);
		//dist.inverseCumulativeProbability(kmSD)
		
		
		// Grid
		kmridRange[0] = kmMean - kmSD*NUMBER_OF_SD_IN_GRID;
		kmridRange[1] = kmMean + kmSD*NUMBER_OF_SD_IN_GRID;
		

		// Density
		double kmdt = (kmridRange[1]-kmridRange[0])/this.ngrid;;
		for (int j = 1; j <= this.ngrid; j ++) {
			
			double kmLower = kmridRange[0] + (j-1) * kmdt;
			double kmUpper = kmridRange[0] + (j) * kmdt;
			
			// Density
			double p1 = this.dist.density(kmLower);
			double p2 = this.dist.density(kmUpper);
			double p = (kmUpper-kmLower) * (p2+p1)/2;
			gridDensity[j] = p;
			
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
		
		double hill = this.getHill();
		if (hill != 1) {
			a = Math.pow(a, hill);
		}
		
		double pSum = 0;
		double vmax = vmaxInput.get().getValue();
		double kmdt = (kmridRange[1]-kmridRange[0])/this.ngrid;;
		
		// Find grid boundaries
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
			double p = gridDensity[j];
			
			
			pSum += p;
			
			// Calculate v
			double vThis = (vmax*a) / (km + a);
			
			//Log.warning("\tp= " + p + " v(a) = v(" + a + ")=" + vThis + " ------ >" + (vThis*p));
			
			v += vThis*p;
			
			
		}
		
		
		v = v / pSum;
		
		//Log.warning("psum = " + pSum + " v(a) = v(" + a + ")=" + v + " | " + kmridRange[0]);
		
		return v;
	}
	

}
