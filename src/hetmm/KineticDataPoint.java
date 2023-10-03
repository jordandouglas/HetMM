package hetmm;


import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;


@Description("A single Michaelis-Menten data pair -a reaction rate (v) and reactant concentration (a)")
public class KineticDataPoint extends BEASTObject {
	
	
	final public Input<String> aNameInput = new Input<>("reactant", "reactant name", Input.Validate.REQUIRED);
	final public Input<Double> aInput = new Input<>("a", "reactant concentration", Input.Validate.REQUIRED);
	final public Input<Double> vInput = new Input<>("v", "reaction velocity concentration", Input.Validate.OPTIONAL);

	
	@Override
	public void initAndValidate() {
		
		if (aInput.get() < 0) {
			throw new IllegalArgumentException("Please ensure a is >= 0");
		}
		
		if (vInput.get() != null && vInput.get() < 0) {
			throw new IllegalArgumentException("Please ensure v is >= 0");
		}
		
	}

	
	public String getReactant() {
		return aNameInput.get();
	}
	
	public double getA() {
		return aInput.get();
	}
	
	public boolean hasV() {
		return vInput.get() != null;
	}
	
	public double getV() {
		return vInput.get();
	}
}
