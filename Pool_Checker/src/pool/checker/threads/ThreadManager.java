/******************************************************************************************************* Copyright (c) 2011 Regents of the University of California.
 * All rights reserved.
 * This software was developed at the University of California, Irvine.
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in
 * the documentation and/or other materials provided with the
 * distribution.
 * 3. All advertising materials mentioning features or use of this
 * software must display the following acknowledgment:
 * "This product includes software developed at the University of
 * California, Irvine by Nicolas Oros, Ph.D.
 * (http://www.cogsci.uci.edu/~noros/)."
 * 4. The name of the University may not be used to endorse or promote
 * products derived from this software without specific prior written
 * permission.
 * 5. Redistributions of any form whatsoever must retain the following
 * acknowledgment:
 * "This product includes software developed at the University of
 * California, Irvine by Nicolas Oros, Ph.D.
 * (http://www.cogsci.uci.edu/~noros/)."
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 * WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 * IN NO EVENT SHALL THE UNIVERSITY OR THE PROGRAM CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. *******************************************************************************************************/
package pool.checker.threads;

import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledExecutorService;

import pool.checker.MainActivity;
import android.os.Handler;
import android.os.Looper;
import android.os.Message;

/** This class simply holds all of the other threads, such as that managing GPS, sensors,
 * or the camera */
public class ThreadManager extends Thread {
	private final String		TAG			= getClass().getSimpleName();
	public Handler mHandler;
	private MainActivity		app;

	private LocationResolver	the_gps;
	private SensorFusion		the_sensors;
	private Solar_Thread		the_sun;

	ScheduledExecutorService	executor	= Executors.newScheduledThreadPool(2);


	public void run() {
        Looper.prepare();

        mHandler = new Handler() {
            public void handleMessage(Message msg) {
                // process incoming messages here
            }
        };

        Looper.loop();
	}
	
	/** Get a new ThreadManager, creating instances of all the threads, and starting them
	 * 
	 * @param app
	 *            - The host application
	 * @param ip
	 *            - The ip of the server to connect to */
	public ThreadManager(MainActivity app) {

		this.app = app;
		the_gps = (LocationResolver) new LocationResolver(app).execute();
		
//		executor.scheduleAtFixedRate(the_gps = new Location_thread(app), 5, 5, TimeUnit.SECONDS);
			
		the_sensors = new SensorFusion(this);
		the_sensors.addResultCallback(app);

//		the_sun = new Solar_Thread(app.getLocation());
//		the_sun.addResultCallback(app);
		
//		the_weather = (ForecastFetcherTask) new ForecastFetcherTask(app, location).execute();
	}

	/** Stop all the threads */
	public void stopAll() {
		// the_cam.stop_thread();
		executor.shutdown();
		// the_sun.stop_thread();
		the_sensors.Stop();
	}

	/** Restart all the threads. */
	public void restartAll() {
	}

	public LocationResolver getGPSThread() {
		return this.the_gps;
	}

	public Solar_Thread getSolarThread() {
		return this.the_sun;
	}

	public SensorFusion getSensorsThread() {
		return this.the_sensors;
	}

	public MainActivity getMainActivity() {
		return this.app;
	}
}
