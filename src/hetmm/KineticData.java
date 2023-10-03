package hetmm;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.IntegerParameter;


@Description("Michaelis-Menten data matrix - reaction rates (v) and reactant concentrations (a)")
public class KineticData extends BEASTObject {

	final public Input<List<KineticDataPoint>> observationsInput = new Input<>("observation", "an observation (v and a)", new ArrayList<>());
	
	
	
	int numObservations = 0;
	
	@Override
	public void initAndValidate() {
		this.numObservations = observationsInput.get().size();
		if (this.numObservations < 2) {
			throw new IllegalArgumentException("Please provide at least 2 observations");
		}
		

		
	}
	
	
	public KineticDataPoint getObservation(int i) {
		return observationsInput.get().get(i);
	}
	
	public int getNumberOfObservations() {
		return this.numObservations;
	}


	public double getA(int i) {
		return this.getObservation(i).getA();
	}
	
	public double getV(int i) {
		return this.getObservation(i).getV();
	}


	public boolean hasV(int i) {
		return this.getObservation(i).hasV();
	}
	

}
