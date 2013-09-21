/**
 * 
 */
package pool.checker.threads;

import java.util.Calendar;
import java.util.Collections;
import java.util.Date;
import java.util.HashSet;
import java.util.Set;
import java.util.Timer;

import org.joda.time.DateTime;
import org.opencv.core.Point;

import pool.checker.MainActivity;
import pool.checker.listeners.MySunListener;
import android.location.Location;
import android.os.Handler;
import android.util.Log;

/** @author Pascal */

public class Solar_Thread implements Runnable {
    /**
     * Nombre de millisecondes dans une journée. Cette constante est
     * utilisée pour convertir des intervalles de temps du Java en
     * nombre de jours.
     */
    private static double MILLIS_IN_DAY = 1000 * 60 * 60 * 24;

    /**
     * Jour julien correspondant à l'époch du Java (1er janvier 1970 à minuit).
     * Cette constante est utilisée pour convertir des dates du Java en jour
     * julien.
     * La valeur {@link #julianDay} du 1er janvier 2000 00:00 GMT est 2451544.5 jours.
     * La valeur {@link Date#getTime} du 1er janvier 2000 00:00 GMT est 10957 jours.
     */
    private static double JULIAN_DAY_1970 = 2451544.5 - 10957;

    private final String TAG = getClass().getSimpleName();

    private Timer solarTimer = new Timer();
    private Handler mSunHandler;

    public static final int TIME_CONSTANT = 30;

    public Point mSunPos;
    private MySunListener mListener;
    private ThreadManager threadManager;
    private DateTime dateTime;

    private Calendar calendar;

    private double min;
    private double sec;
    private double heu;

    private Thread processingThread;
    private MainActivity app;

    private boolean STOP_THREAD;

    public Location location;

    public boolean loc_flag = false;

    private Set<ResultCallback> mResultCallbacks = Collections
	    .synchronizedSet(new HashSet<ResultCallback>());

    public Solar_Thread(Location location) {
	// this.threadManager = location;
	this.location = location;
	mSunHandler = new Handler();
	// solarTimer.scheduleAtFixedRate(new calculateSunTask(), (long)1000, (long)TIME_CONSTANT);

	// init();
    }

    public Solar_Thread(ThreadManager threadManager) {
	// TODO Auto-generated constructor stub
	this.location = threadManager.getGPSThread().location;
    }

    public Solar_Thread(MainActivity app) {
	// TODO Auto-generated constructor stub
	// TODO Auto-generated constructor stub
	// this.location = threadManager.getGPSThread().location;
    }

    public Point soleil(Location location/* $nebu */) {

	/* $cloud = $nebu / 1024;
	 * // on va dire que c'est l’indice de Perreaudau (In).
	 * // Pour les conditions claires, In>0.8. Pour les conditions très maussades,
	 * In<0.1.
	 * // re.jrc.ec.europa.eu/pvgis/solres/solmod3.htm
	 * $clientrawfile = 'modules/meteo/clientraw.txt';
	 * $dataraw = file_get_contents($clientrawfile);
	 * // clean up any blank lines
	 * $dataraw = trim($dataraw);
	 * $dataraw = preg_replace("/[\r\n]+[\s\t]*[\r\n]+/","\n",$dataraw);
	 * $data = explode(" ", $dataraw);
	 * // TEMPERATURE
	 * (float)$TempC = $data[4];
	 * // HUMIDITY
	 * (int)$Humidity = $data[5];
	 * // WIND
	 * $Wind = number_format($data[1] / 1.852, 1, '.', ''); // kts
	 * $Gust = number_format($data[2] / 1.852, 1, '.', ''); // kts
	 * (int)$WindDir = $data[3];
	 * // BAROMETER
	 * (float)$Baro = $data[6];
	 * // RAIN
	 * (float)$RainToday = $data[7]; */

	double LAT = deg2rad(location.getLatitude()); // deg2rad(42.7136788229883);
	double LON = deg2rad(location.getLongitude()); // deg2rad(2.84930393099784);

	double SO = 1367.6;
	// SET @JOU = DAYOFYEAR(NOW());
	calendar = Calendar.getInstance();
	Date now = new Date();
	dateTime = new DateTime(now);
	heu = calendar.get(Calendar.HOUR_OF_DAY);
	min = calendar.get(Calendar.MINUTE);
	sec = calendar.get(Calendar.SECOND);

	// double DST = time.gmtoff /3600;
	double DST = (calendar.get(Calendar.ZONE_OFFSET) + calendar
		.get(Calendar.DST_OFFSET)) / 60 / 60 / 1000;
	/* double JOU = calendar.get(Calendar.DAY_OF_YEAR); // date("z");
	 * // double HOD = (date("G")*60 + date("i") - .5) / 60; //G Heure format 24 i
	 * // minutes
	 * double HOD = (calendar.get(Calendar.HOUR_OF_DAY) * 60
	 * + calendar.get(Calendar.MINUTE) - .5) / 60;
	 * double HRA = 2 * pi() * (HOD - 12.0) / 24.0; // Angle Horaire
	 * double jd = GregorianToJD(calendar.get(Calendar.YEAR),
	 * calendar.get(Calendar.MONTH), calendar.get(Calendar.DAY_OF_MONTH));
	 * // double HD = ((date("G")+((date('i') + date('s') / 60)/60))/24);
	 * double HD = ((calendar.get(Calendar.HOUR_OF_DAY) + ((calendar
	 * .get(Calendar.MINUTE) + calendar.get(Calendar.SECOND) / 60) / 60)) / 24);
	 * // correct for half-day offset
	 * // double dayfrac = date('G') / 24 - .5;
	 * double dayfrac = calendar.get(Calendar.HOUR_OF_DAY) / 24 - .5;
	 * if (dayfrac < 0)
	 * dayfrac += 1;
	 * // now set the fraction of a day
	 * // double frac = dayfrac + (date('i') + date('s') / 60) / 60 / 24;
	 * double frac = dayfrac
	 * + (calendar.get(Calendar.MINUTE) + calendar.get(Calendar.SECOND) / 60)
	 * / 60 / 24;
	 * double julianDate = jd + frac - DST / 24;
	 * // echo "JD $julianDate<br>";
	 * julianDate = julianDay(new Date().getTime()) - DST /24;
	 * //julianCentury = julianCentury(new Date());
	 * double julianCentury = (julianDate - 2451545) / 36525; // G2
	 * // echo "JC $julianCentury<br>"; */
	double julianCentury = julianCentury(now);

	// Log.e(TAG, "JulianCentury " + julianCentury + " DST " + DST);

	double GMLS = fmod(280.46646 + julianCentury
		* (36000.76983 + julianCentury * 0.0003032), 360); // Geom Mean Long Sun
								   // (deg) I2
	double GMAS = 357.52911 + julianCentury
		* (35999.05029 - 0.0001537 * julianCentury); // Geom Mean Anom Sun (deg)

	// J2
	// echo "GMLS $GMLS/GMAS $GMAS<br>";
	double EEO = 0.016708634 - julianCentury
		* (0.000042037 + 0.0000001267 * julianCentury); // Eccent Earth Orbit K2
	double SEC = SIN(deg2rad(GMAS))
		* (1.914602 - julianCentury
			* (0.004817 + 0.000014 * julianCentury))
		+ SIN(deg2rad(2 * GMAS))
		* (0.019993 - 0.000101 * julianCentury)
		+ SIN(deg2rad(3 * GMAS)) * 0.000289;// Sun Eq of Ctr L2
	double STL = GMLS + SEC; // Sun True Long (deg) M2
	double STA = GMAS + SEC; // Sun True Anom (deg) N2
	double SRV = (1.000001018 * (1 - EEO * EEO))
		/ (1 + EEO * COS(deg2rad(STA))); // Sun
						 // Rad
						 // Vector
						 // (AUs)
						 // O2
	// echo "SRV $SRV<br>";
	double SAL = STL - 0.00569 - 0.00478
		* SIN(deg2rad(125.04 - 1934.136 * julianCentury)); // Sun App Long (deg)
								   // P2
	double MOE = 23 + (26 + ((21.448 - julianCentury
		* (46.815 + julianCentury
			* (0.00059 - julianCentury * 0.001813)))) / 60) / 60; // Mean
									      // Obliq
									      // Ecliptic
									      // (deg)
									      // Q2
	double OC = MOE + 0.00256
		* COS(deg2rad(125.04 - 1934.136 * julianCentury)); // Obliq
								   // Corr
								   // (deg)
								   // R2
	double SRA = rad2deg(ATAN2(COS(deg2rad(OC)) * SIN(deg2rad(SAL)),
		COS(deg2rad(SAL)))); // Sun Rt Ascen (deg) S2
	double SD = rad2deg(ASIN(SIN(deg2rad(OC)) * SIN(deg2rad(SAL)))); // Sun Declin
									 // (deg) T2 ->
									 // DES
	// echo "SD $SD<br>";
	double VY = TAN(deg2rad(OC / 2)) * TAN(deg2rad(OC / 2)); // var y U2
	// echo "SRA $SRA VY $VY<br>";
	// Eq of Time (minutes) V2
	double EQT = 4 * rad2deg(VY * SIN(2 * deg2rad(GMLS)) - 2 * EEO
		* SIN(deg2rad(GMAS)) + 4 * EEO * VY * SIN(deg2rad(GMAS))
		* COS(2 * deg2rad(GMLS)) - 0.5 * VY * VY
		* SIN(4 * deg2rad(GMLS)) - 1.25 * EEO * EEO
		* SIN(2 * deg2rad(GMAS)));
	// echo "EQT $EQT<br>";
	// Log.e(TAG, "SD" + SD + " VY " + VY + " EQT "+ EQT); //OK

	
	// DEGRES(ACOS(COS(RADIANS(90,833))/(COS(RADIANS($B$3))*COS(RADIANS(T136)))-TAN(RADIANS($B$3))*TAN(RADIANS(T136))))
	/** @ HA Sunrise (deg) W2 */
	double HASR = rad2deg(ACOS(COS(deg2rad(90.833))
		/ (COS((LAT)) * COS(deg2rad(SD))) - TAN((LAT))
		* TAN(deg2rad(SD))));
	// echo "HASR $HASR<br>";
	double QN = (720 - 4 * rad2deg(LON) - EQT + DST * 60) / 1440;// Solar Noon
								     // (LST)
								     // X2
	double SRT = QN - HASR * 4 / 1440; // Sunrise Time (LST) Y2
	double SST = QN + HASR * 4 / 1440; // Sunset Time (LST) Z2
	double SLD = 8 * HASR; // Sunlight Duration (minutes) AA2

	// double HD = ((date("G")+((date('i') + date('s') / 60)/60))/24);
	// =(15+(35+50/60)/60)/24
	// Log.w(TAG, heu + " " + min + " " + sec);
	double HD = (heu + (min + sec / 60) / 60) / 24;
	double TST = fmod(HD * 1440 + EQT + 4 * rad2deg(LON) - 60 * DST, 1440); // True
										// Solar
										// Time
										// (min)
										// AB2

	double HA = (TST / 4 < 0) ? TST / 4 + 180 : TST / 4 - 180; // Hour Angle (deg) AC2
	double SZA = rad2deg(ACOS(SIN((LAT)) * SIN(deg2rad(SD)) + COS((LAT))
		* COS(deg2rad(SD)) * COS(deg2rad(HA))));// Solar_Thread Zenith Angle (deg)
	// echo "TST $TST / SZA $SZA<br>";
	double SEA = 90 - SZA; // Solar_Thread Elevation Angle (deg) AE2
	// Log.e(TAG, "HD " + HD + " HA" + HA + " SZA " + SZA + " SEA " + SEA);

	// Approx Atmospheric Refraction (deg) AF2
	// switch(SEA) {
	double AAR = 0;
	if (SEA > 85 && SEA < 90)
	    AAR = 0;
	else if (SEA < -0.575)
	    AAR = (-20.772 / TAN(deg2rad(SEA))) / 3600;
	else if (SEA > -0.575 && SEA < 5)
	    AAR = (1735 + SEA
		    * (-518.2 + SEA * (103.4 + SEA * (-12.79 + SEA * 0.711)))) / 3600;
	else if (SEA > 5 && SEA < 85)
	    AAR = (58.1 / TAN(deg2rad(SEA)) - 0.07
		    / Math.pow(TAN(deg2rad(SEA)), 3) + 0.000086 / Math.pow(
		    TAN(deg2rad(SEA)), 5)) / 3600;
	// }
	// echo "SEA: $SEA / AAR: $AAR<br>";
	SEC = SEA + AAR; // Solar_Thread Elevation corrected for atm refraction (deg)
	// Solar_Thread Azimuth Angle (deg cw from N)
	double SAA = (HA > 0) ? fmod(
		rad2deg(ACOS(((SIN((LAT)) * COS(deg2rad(SZA))) - SIN(deg2rad(SD)))
			/ (COS((LAT)) * SIN(deg2rad(SZA))))) + 180, 360)
		: fmod(540 - rad2deg(ACOS(((SIN((LAT)) * COS(deg2rad(SZA))) - SIN(deg2rad(SD)))
			/ (COS((LAT)) * SIN(deg2rad(SZA))))), 360);
	// echo "SEC: $SEC / SAA: $SAA<br>";

	/* FE = 1. + 0.0334*COS(2.*PI()*(JOU-2.7206)/365.25);
	 * //$DES = 0.4093*SIN($JOU*((2.*PI())/365.)-1.405);
	 * //$H = (sin($LAT)*sin($DES)+cos($LAT)*cos($DES)*cos($HRA));
	 * //$HS = rad2deg($H);
	 * H = deg2rad(SEC);
	 * SOL = SO*FE*H;
	 * //pow( $a , 0.6 ) use: exp( 0.6 * log($a) )
	 * Alt = 54;
	 * //$PAtm = (1013.25 * (1 - (0.0065 * $Alt) / 288.15)^5.255) / 100; // en Pa
	 * Pvs = 2.165 * exp(8.02 * log(1.098 + TempC / 100)); // en mmHg (millimètre de mercure)
	 * Pv = Pvs * Humidity / 100;
	 * //$m = $PAtm / (101325 * sin($HS) + 15198.75 * (3.885 - $HS)^(-1.253)); //masse d'air
	 * optique relative (m) h est la hauteur du soleil en degrés
	 * //$m = (0.89^$Alt) / sin($H);
	 * // http://www.cder.dz/download/Art9-4_5.pdf
	 * // $m = (1 - $Alt/10000) / (sin($H) + 0.50572 * exp(-1.6364 * log(6.07995 + $H))); //
	 * formule 3
	 * // http://www.cythelia.fr/images/file/Gisement-solaire_Alain%20Ricaud_Jan-2011.pdf
	 * m = exp((-Alt/1000)/8434.5) / (sin(H) + 0.50572 * exp(-1.6364 * log(6.07995 + H))); //
	 * formule 3
	 * //echo "MA".$m." ".(1/cos(deg2rad($SZA)));
	 * //$ER = 1 / (0.9 * $m + 9.4); //épaisseur optique de Rayleigh (ER)
	 * // http://www.cder.dz/download/Art9-4_5.pdf
	 * ER = 1/(6.55567 + 1.7513 * m - 0.1202 * pow(m,2) + 0.0065 * pow(m,3) - 0.00013 *
	 * pow(m,4));
	 * //echo "RA".$ER." ";
	 * B = 0.075; // pour un lieu rural 0,02 pour un lieu situé en montagne 0,20 pour un lieu
	 * industriel (atmosphère polluée)
	 * TL = 2.4 + 14.6 * B + 0.4 * (1 + 2 * B) * log(Pv); //facteur de trouble de Linke
	 * SOL = SOL * exp(-ER * $m * TL); // rayonnement solaire direct sur un plan récepteur
	 * normal à ce rayonnement en W/m²
	 * if (SOL <= 0 || $SOL > 3000) return "0.0";
	 * //Albédo 0.22 //0 pour un plan horizontal 90 pour vertical
	 * DI =
	 * 125*exp(0.4*log(sin(H)))*((1+cos(deg2rad(0)))/2)+211.86*exp(1.22*log(sin(H)))*((1-cos
	 * (deg2rad(0)))/2);
	 * if(is_nan(DI)) DI=0;
	 * //echo "DI".$DI;
	 * CI = 1; // Coefficient d'orientation à calculer
	 * SOL = (cloud * SOL*CI) + DI /* /* $cloud; */
	mSunPos = new Point(SAA, SEC);
	// update sensor output in GUI
	// mSunHandler.post(updateSunOrientation);
	// notifyResultCallback(mSunPos);
	// Log.e(TAG, "Solar" + mSunPos.toString());
	return mSunPos;
	// http://herve.silve.pagesperso-orange.fr/solaire.htm
	// $DSO = ACOS(-TAN(($LAT))*TAN($DES));
	// return number_format(SOL, 1, '.');
	// SET @LAT = RADIANS(42.713674); #A prendre dans la table station
	// SET @S0 = 1367.6; #Wm/2 ou 15.392 pour mm/jour
	// SET @HOD = (HOUR(NOW())*60+MINUTE(NOW())-.5)/60; #Heure decimale
	// SET @HRA = 2*PI()*(@HOD-12.0)/24.0; #the hour angle
	// #Reste Equation du temps pour etre un peu plus juste..

	// SET @FE = 1. + 0.0334*COS(2.*PI()*(@JOU-2.7206)/365.25); #the excentricity
	// factor or the relative distance between Earth and Sun VERIFIE
	// SET @DES = 0.4093*SIN(@JOU*((2.*PI())/365.)-1.405); #Delta the Solar_Thread
	// declination in radians VERIFIE

	// SET @SOL = @S0*@FE*(sin(@LAT)*sin(@DES)+cos(@LAT)*cos(@DES)*cos(@HRA)); #

	// #SET @DSO = ACOS(-TAN(RADIANS(@LAT))*TAN(@DES));
	// SET NEW.NoJour = @JOU;
	// SET NEW.solar = IF(@SOL<=0,0,@SOL);
    }

    /**
     * Retourne le jour julien d'une date. Il ne s'agit pas du jour julien dans
     * une année. Ce jour julien lÃ (nommÃ© ainsi pour <i>Julius Scaliger</i>, et
     * non <i>Julius Caesar</i>) est le nombre de jours Ã©coulÃ©s depuis midi le
     * 1er janvier 4713 avant JÃ©sus-Christ.
     */
    public static double julianDay(final Date time) {
	return julianDay(time.getTime());
    }

    /**
     * Computes the {@linkplain #julianDay(Date) julian day}.
     * 
     * @param time
     *        The date in milliseconds ellapsed since January 1st, 1970.
     */
    static double julianDay(final long time) {
	return (time / MILLIS_IN_DAY) + JULIAN_DAY_1970;
    }

    /**
     * Retourne le nombre de siÃ¨cles Ã©coulÃ©s depuis le 1 janvier 2000 Ã midi.
     * Cette information est utilisÃ©e dans les formules de Laskar (1986) pour
     * calculer la longueur d'une annÃ©e tropicale, ainsi que par Chapront-Touze
     * et Chapront (1988) pour la longueur d'un mois synodique.
     */
    static double julianCentury(final Date time) {
	return ((time.getTime() / MILLIS_IN_DAY) + (JULIAN_DAY_1970 - 2451545.0)) / 36525;
    }

    private double floor(double d) {
	return Math.floor(d);
    }

    private double pi() {
	return Math.PI;
    }

    private double fmod(double d, int i) {
	return d % i;
    }

    private double ACOS(double d) {
	return Math.acos(d);
    }

    private double TAN(double d) {
	return Math.tan(d);
    }

    private double rad2deg(double d) {
	return Math.toDegrees(d);
    }

    private double ASIN(double d) {
	return Math.asin(d);
    }

    private double ATAN2(double x, double y) {
	return Math.atan2(y, x);
    }

    private double COS(double d) {
	return Math.cos(d);
    }

    private double SIN(double d) {
	return Math.sin(d);
    }

    private double deg2rad(double d) {
	return Math.toRadians(d);
    }

    public void stop_thread() {
	STOP_THREAD = true;
    }

    public Solar_Thread start_thread() {
	STOP_THREAD = false;
	processingThread = new Thread(this);
	processingThread.start();
	return null;
    }

    public Point setSunPos(Location location) {
	mSunPos = soleil(location);
	return mSunPos;
    }

    public Point getSunPos() {
	return mSunPos;
    }

    public void updateSun() {
	setSunPos(location);
	Log.i(TAG, mSunPos.toString());
	notifyResultCallback(mSunPos);
    }

    private void notifyResultCallback(Point mSunPos) {
	for (ResultCallback resultCallback : mResultCallbacks) {
	    try {
		resultCallback.onSolarReady(mSunPos);
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
	void onSolarReady(Point mSunPos);
    }

    @Override
    public void run() {
	updateSun();
    }

}
