package pool.checker.threads;

import java.util.List;

import pool.checker.MainActivity;
import android.app.Activity;
import android.app.ProgressDialog;
import android.content.Context;
import android.location.Location;
import android.location.LocationListener;
import android.location.LocationManager;
import android.os.AsyncTask;
import android.os.Bundle;
import android.util.Log;

import com.hourlyweather.forecast.ForecastFetcherTask;

/**
 * Asynchronous task to determine the users location from any available source
 * and inform the user on the progress
 * 
 * @author dhgonsalves
 * 
 */
public class LocationResolver extends AsyncTask<Void, String, Location> {
    private final String TAG = getClass().getSimpleName();
    private static Solar_Thread the_sun;

    private ThreadManager threadManager;
    private MainActivity app;

    protected ProgressDialog dialog;
    // protected ForecastFetcherInt application;
    protected LocationManager locationManager;
    protected Location location;

    public static ForecastFetcherTask the_weather;

    public LocationResolver(Activity application) {
	location = null;
	dialog = new ProgressDialog(application);

	locationManager = (LocationManager) application
		.getSystemService(Context.LOCATION_SERVICE);

	this.app = (MainActivity) application;
    }

    /*public LocationResolver(HourlyWeather application) {
     * this((Activity) application);
     * this.application = application;
     * }
     * 
     * public LocationResolver(HourlyWidgetConfigure application) {
     * this((Activity) application);
     * this.application = application;
     * } */
    protected void onPreExecute() {

	dialog.setMessage("Determining Location");
	dialog.show();

	List<String> providers = locationManager.getAllProviders();
	// set up listening for a location update
	for (String provider : providers)
	    if (locationManager.isProviderEnabled(provider))
		locationManager.requestLocationUpdates(provider, 0, 0,
			locationListener);
    }

    @Override
    protected Location doInBackground(final Void... unused) {

	// check every half second for updates for 10 seconds
	for (int i = 0; true; i++) {
	    if (location == null)
		try {
		    Thread.sleep(500);
		} catch (InterruptedException e) {
		    // can't do anything here but continue
		}
	    else
		break;

	    // if the user has waited more then 2 seconds tell them we are still
	    // waiting
	    if (i > 5)
		publishProgress("waiting for location.");
	    else if (i > 10)
		publishProgress("still waiting for location.");
	}

	return location;
    }

    LocationListener locationListener = new LocationListener() {
	public void onLocationChanged(Location currentLocation) {
	    location = currentLocation;
	    locationManager.removeUpdates(this);
	}

	public void onProviderDisabled(String provider) {
	}

	public void onProviderEnabled(String provider) {
	}

	public void onStatusChanged(String provider, int status, Bundle extras) {
	}
    };

    @Override
    protected void onProgressUpdate(String... values) {
	dialog.setMessage(values[0]);
    }

    @Override
    protected void onPostExecute(Location location) {
	locationManager.removeUpdates(locationListener);

	if (location == null) {
	    // notify user that their location isn't availiable
	    dialog.setMessage("Your location is currently unavailable, try again later.");
	    return;
	}

	// fetch the forecast for the location
	Log.w(TAG, location.toString());
	// app.location.set(location);

	the_sun = new Solar_Thread(location);
	the_sun.addResultCallback(app);
	the_sun.start_thread();
	if (deviceIsReadyToPullForecast(this, locationManager)) {
	    the_weather = (ForecastFetcherTask) new ForecastFetcherTask(app,
		    location).execute();
	    the_weather.addResultCallback(app);
	}
	dialog.dismiss();

    }

    private boolean deviceIsReadyToPullForecast(
	    LocationResolver locationResolver, LocationManager locationManager) {
	// TODO Auto-generated method stub
	if (!com.hourlyweather.yrno.NetworkUtil.isNetworkAvailable(app))
	    return false;
	return true;
    }
}
