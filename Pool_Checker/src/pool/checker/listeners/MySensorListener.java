package pool.checker.listeners;

public interface MySensorListener {

	public void gotNewAccelerometerData(float x, float y, float z);
	public void gotNewOrientationData(float x, float y, float z);
	
	public void gotNewDisplacement(float cap, float dep, float vel);
}
