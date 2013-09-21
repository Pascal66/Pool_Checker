package com.hourlyweather.forecast;

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

import android.app.ProgressDialog;
import android.content.Context;
import android.location.Location;
import android.os.AsyncTask;

/**
 * task that handles asynchronous pulling of a forecast as well as notifying
 * user of progress
 * 
 * @author dhgonsalves
 */
public class ForecastFetcherTask extends
	AsyncTask<Void, String, HourlyForecast> {

    private final String TAG = getClass().getSimpleName();

    protected Context context;
    protected ProgressDialog loadingDialog;
    protected Location location;

    private Set<ResultCallback> mResultCallbacks = Collections
	    .synchronizedSet(new HashSet<ResultCallback>());

    public ForecastFetcherTask(Context context) {
	this.context = context;
	loadingDialog = new ProgressDialog(context);
	loadingDialog.setCancelable(true);

	location = null;
    }

    public ForecastFetcherTask(Context context, Location location) {
	this(context);
	this.location = location;
    }

    @Override
    protected void onPreExecute() {
	loadingDialog.setMessage("Loading Forecast");
	loadingDialog.show();
    }

    @Override
    protected HourlyForecast doInBackground(Void... v) {

	// publishProgress("Loading Forecast");
	HourlyForecast forecast = ForecastCacher.getForecast(context, location);

	if (forecast != null)
	    publishProgress("Prepairing Forecast");

	loadingDialog.dismiss();
	// if there were a network related error forecast will be null
	return forecast;
    }

    @Override
    protected void onPostExecute(HourlyForecast forecast) {
	notifyResultCallback(forecast);
	loadingDialog.dismiss();
    }

    private void notifyResultCallback(HourlyForecast forecast) {
	for (ResultCallback resultCallback : mResultCallbacks) {
	    try {
		resultCallback.onWeatherReady(forecast);
	    } catch (Exception e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	    }
	}
    }

    public void addResultCallback(ResultCallback resultCallback) {
	mResultCallbacks.add(resultCallback);
    }

    public void removeResultCallback(ResultCallback resultCallback) {
	mResultCallbacks.remove(resultCallback);
    }

    public interface ResultCallback {
	void onWeatherReady(HourlyForecast forecast);
    }
}
