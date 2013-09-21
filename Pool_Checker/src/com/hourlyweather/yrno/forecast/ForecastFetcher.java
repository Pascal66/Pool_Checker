package com.hourlyweather.yrno.forecast;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.StringReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.text.ParseException;

import org.joda.time.DateTime;
import org.joda.time.DateTimeZone;
import org.joda.time.Hours;
import org.joda.time.MutableDateTime;
import org.joda.time.format.DateTimeFormat;
import org.joda.time.format.DateTimeFormatter;
import org.xmlpull.v1.XmlPullParser;
import org.xmlpull.v1.XmlPullParserException;
import org.xmlpull.v1.XmlPullParserFactory;

import android.util.Log;

import com.hourlyweather.forecast.ForecastHour;
import com.hourlyweather.forecast.HourlyForecast;
import com.hourlyweather.yrno.XmlParserUtil;

public class ForecastFetcher {
    private static final DateTimeFormatter df = DateTimeFormat.forPattern(
	    "YYYY-MM-dd'T'HH':00:00Z'").withZone(DateTimeZone.UTC);

    @SuppressWarnings("unused")
    private final String TAG = getClass().getSimpleName();

    /**
     * Connects to the yr.no api and returns an input stream for the hourly
     * forecast xml
     * 
     * @param lat the latitude your polling for
     * @param lon the longitude your polling for
     * @return
     */
    private static InputStream getHourlyForecastData(HourlyForecast forecast) {

	// pull the weather xml
	try {
	    URL weatherUrl = new URL(
		    "http://api.yr.no/weatherapi/locationforecast/1.8/?lat="
			    + forecast.getLat() + ";lon=" + forecast.getLon());
	    // Log.i(" URL :", weatherUrl+"");
	    return weatherUrl.openStream();
	} catch (MalformedURLException e) {
	    // TODO Auto-generated catch block
	    e.printStackTrace();
	} catch (IOException e) {
	    // TODO Auto-generated catch block
	    e.printStackTrace();
	}

	return null;
    }
    /**
     * Connects to the noaa api and returns an input stream neares airport xml
     * @param lat the latitude your polling for
     * @param lon the longitude your polling for
     * @return
     */
    private static InputStream getNearestAirportData(HourlyForecast forecast) {

	// pull the result from Noaa xml
	try {
	    URL AirportUrl = new URL(
		    "http://www.ncdc.noaa.gov/cdo-services/services/datasets/GHCND/locationsearch?" +
		    		"latitude="+ forecast.getLat()+"&longitude="+ forecast.getLon() +
		    		"&token=kClxOskPZexpKAcoLJeipIAUvGBImaMI&radius=50");
	    // Log.i(" URL :", weatherUrl+"");
	    return AirportUrl.openStream();
	} catch (MalformedURLException e) {
	    // TODO Auto-generated catch block
	    e.printStackTrace();
	} catch (IOException e) {
	    // TODO Auto-generated catch block
	    e.printStackTrace();
	}
	return null;
    }
    /**
     * gets the current hourly forecast based on the supplied latitude and
     * longitude
     * 
     * @param lat
     * @param lon
     * @param hours the number of hours into the future to poll for
     * @return
     */
    public static boolean getHourlyForecast(HourlyForecast forecast) {
	// add the sun rise/set times to the forecast
	// SunriseFetcher.addSunlightDurationToForecast(forecast);

	InputStream input = null; 
	try {
	    input = getHourlyForecastData(forecast);
	} catch (Exception e) {
	    // TODO Auto-generated catch block
	    e.printStackTrace();
	    return false;
	}
	
	DateTime endWindow;
	{
	    // add a day plus our start time to create a cut off time
	    MutableDateTime temp = forecast.getStart().toMutableDateTime();
	    temp.addHours(forecast.getHours());
	    temp.addDays(1);
	    endWindow = temp.toDateTime();
	}

	XmlPullParser xpp;
	try {
	    XmlPullParserFactory factory = XmlPullParserFactory.newInstance();
	    factory.setNamespaceAware(true);
	    xpp = factory.newPullParser();
	} catch (XmlPullParserException e) {
	    System.out.println("error creating xml parser: " + e.getMessage());
	    return false;
	}

	try {
	    xpp.setInput(new BufferedReader(new InputStreamReader(input)));
	    int eventType = xpp.getEventType();
	    while (eventType != XmlPullParser.END_DOCUMENT) {
		if (eventType == XmlPullParser.START_TAG)
		    if ("time".equalsIgnoreCase(xpp.getName()))
			if (!parseTimeTag(xpp, endWindow, forecast))
			    break;

		eventType = xpp.next();

	    }
	    input.close();
	} catch (Exception e) {
	    // TODO Auto-generated catch block
	    e.printStackTrace();
	    return false;
	}

	return true;
    }
    /**
     * gets the nearest airport based on the supplied latitude and
     * longitude
     * 
     * @param lat
     * @param lon
     * @param hours the number of hours into the future to poll for
     * @return
     */
    public static boolean getNearestAirport(HourlyForecast forecast) {
	InputStream input = null;
	try {
	    input = getNearestAirportData(forecast);
	} catch (Exception e) {
	    // TODO Auto-generated catch block
	    e.printStackTrace();
	    return false;
	}

	DateTime endWindow;
	{
	    // add a day plus our start time to create a cut off time
	    MutableDateTime temp = forecast.getStart().toMutableDateTime();
	    temp.addHours(forecast.getHours());
	    temp.addDays(1);
	    endWindow = temp.toDateTime();
	}

	XmlPullParser xpp;
	try {
	    XmlPullParserFactory factory = XmlPullParserFactory.newInstance();
	    //factory.setNamespaceAware(true);
	    xpp = factory.newPullParser();
	} catch (XmlPullParserException e) {
	    System.out.println("error creating xml parser: " + e.getMessage());
	    return false;
	}

	try {
	    xpp.setInput(new BufferedReader(new InputStreamReader(input)));
	    
	    int eventType = xpp.getEventType();
	    
	    while (eventType != XmlPullParser.END_DOCUMENT) {	
		Log.e("StringReader",xpp.getEventType()+" ");	
		if (eventType == XmlPullParser.START_TAG) 
		    
		    //if ("searchResult".equalsIgnoreCase(xpp.getName()))
			//if (!parseAirportTag(xpp, forecast))
			    //break;

		eventType = xpp.next();

	    }
	    input.close();
	} catch (Exception e) {
	    // TODO Auto-generated catch block
	    e.printStackTrace();
	    return false;
	}

	return true;
    }
/*
<searchResultCollection pageCount="1" totalCount="1">
<searchResult>
<id>GHCND:FR000007747</id>
<minDate>1901-01-01</minDate>
<maxDate>2013-08-19</maxDate>
<name>PERPIGNAN, FR</name>
<score>40</score>
<type>station</type>
<number>1</number>
<inCart>false</inCart>
<inDateRange>false</inDateRange>
</searchResult>
</searchResultCollection>
*/
    private static boolean parseAirportTag(XmlPullParser xpp, HourlyForecast forecast) throws XmlPullParserException, IOException {
		// parse the child elements of the searchResult element for the airport details
		ForecastHour forecastHour = new ForecastHour();
		for (int eventType = xpp.next(); !(eventType == XmlPullParser.END_TAG && "searchResult"
			.equalsIgnoreCase(xpp.getName())); eventType = xpp.next())

		    if (eventType == XmlPullParser.START_TAG)
Log.i("Airport xml", xpp.getName());
			// find all the airport data
			if ("id".equalsIgnoreCase(xpp.getName())) {
			    forecastHour.setAirportId(xpp.getText()/*XmlParserUtil.getAttributeByName(
				    xpp, "id")*/);
			} else if ("name".equalsIgnoreCase(xpp.getName())) {
			    forecastHour.setAirportName(xpp.getText()/*XmlParserUtil.getAttributeByName(xpp,
				    "value")*/);
			}

		//forecast.add(from, to, forecastHour);

		return true;
    }
    /**
     * returns if the current XML pull parser tag is of type "time" and the from
     * attribute is outside of our range
     * 
     * @param xpp
     * @param startHour
     * @param endHour
     * @return
     * @throws IOException
     * @throws XmlPullParserException
     * @throws ParseException
     */
    private static boolean parseTimeTag(XmlPullParser xpp, DateTime endWindow,
	    HourlyForecast forecast) throws XmlPullParserException,
	    IOException, ParseException {

	DateTime from = df.parseDateTime(
		XmlParserUtil.getAttributeByName(xpp, "from")).withZone(
		forecast.getTimeZone());

	// check if we are at the end of our time range
	if (from.isAfter(endWindow))
	    return false;

	DateTime to = df.parseDateTime(
		XmlParserUtil.getAttributeByName(xpp, "to")).withZone(
		forecast.getTimeZone());

	// parse the child elements of the time element for the forecast details
	ForecastHour forecastHour = new ForecastHour();
	for (int eventType = xpp.next(); !(eventType == XmlPullParser.END_TAG && "time"
		.equalsIgnoreCase(xpp.getName())); eventType = xpp.next())

	    if (eventType == XmlPullParser.START_TAG)

		// find all the weather symbols
		if ("symbol".equalsIgnoreCase(xpp.getName()))
		    forecastHour.setSymbol(XmlParserUtil.getAttributeByName(
			    xpp, "id"));
		else if ("precipitation".equalsIgnoreCase(xpp.getName())) {
		    Double precipitation = getPrecipitation(from, to,
			    XmlParserUtil.getAttributeByName(xpp, "value"));
		    forecastHour.setPrecipitation(precipitation);
		} else if ("windspeed".equalsIgnoreCase(xpp.getName())) {
		    Double windSpeed = getWindSpeed(from, to,
			    XmlParserUtil.getAttributeByName(xpp, "mps"));
		    forecastHour.setWindSpeed(windSpeed);
		} else if ("temperature".equalsIgnoreCase(xpp.getName())) {
		    forecastHour.setTemp(XmlParserUtil.getAttributeByName(xpp,
			    "value"));
		}

	forecast.add(from, to, forecastHour);

	return true;
    }

    private static Double getWindSpeed(DateTime from, DateTime to,
	    String windSpeedString) {

	Double windSpeed;
	try {
	    if (windSpeedString != null)
		windSpeed = Double.parseDouble(windSpeedString);
	    else
		return null;
	} catch (NumberFormatException e) {
	    // nothing we can do now, the api may have changed
	    return null;
	}

	return windSpeed;
    }

    /**
     * gets the percipitation per hour based on the span and the converted
     * precipitation string
     * 
     * @param from
     * @param to
     * @param precipitationString
     * @return
     */
    private static Double getPrecipitation(DateTime from, DateTime to,
	    String precipitationString) {
	Double precipitation;
	try {
	    if (precipitationString != null)
		precipitation = Double.parseDouble(precipitationString);
	    else
		return null;
	} catch (NumberFormatException e) {
	    // nothing we can do now, the api may have changed
	    return null;
	}

	return precipitation / Hours.hoursBetween(from, to).getHours();
    }
}
