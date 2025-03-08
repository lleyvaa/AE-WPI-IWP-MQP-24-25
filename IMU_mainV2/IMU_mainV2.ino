//Lyle Edwards and Zack Stahl
//Feb_04_2025
/*
The sensor outputs provided by the library are the raw 16-bit values
obtained by concatenating the 8-bit high and low accelerometer and
magnetometer data registers. They can be converted to units of g and
gauss using the conversion factors specified in the datasheet for your
particular device and full scale setting (gain).

Example: An LSM303D gives a magnetometer X axis reading of 1982 with
its default full scale setting of +/- 4 gauss. The M_GN specification
in the LSM303D datasheet (page 10) states a conversion factor of 0.160
mgauss/LSB (least significant bit) at this FS setting, so the raw
reading of -1982 corresponds to 1982 * 0.160 = 317.1 mgauss =
0.3171 gauss.

In the LSM303DLHC, LSM303DLM, and LSM303DLH, the acceleration data
registers actually contain a left-aligned 12-bit number, so the lowest
4 bits are always 0, and the values should be shifted right by 4 bits
(divided by 16) to be consistent with the conversion factors specified
in the datasheets.

Example: An LSM303DLH gives an accelerometer Z axis reading of -16144
with its default full scale setting of +/- 2 g. Dropping the lowest 4
bits gives a 12-bit raw value of -1009. The LA_So specification in the
LSM303DLH datasheet (page 11) states a conversion factor of 1 mg/digit
at this FS setting, so the value of -1009 corresponds to -1009 * 1 =
1009 mg = 1.009 g.
*/

#include <SD.h>
#include <Adafruit_Sensor.h>
#include <Adafruit_L3GD20_U.h>
#include <Adafruit_LSM303_U.h>
 
File myFile;

char report[200];
Adafruit_L3GD20_Unified gyro = Adafruit_L3GD20_Unified(20);
Adafruit_LSM303_Accel_Unified accel = Adafruit_LSM303_Accel_Unified(30);

void displaySensorDetails(void)
{
  //Gyro
  sensor_t sensorGyro;
  gyro.getSensor(&sensorGyro);
  Serial.println("------------------------------------");
  Serial.print  ("Sensor:       "); Serial.println(sensorGyro.name);
  Serial.print  ("Driver Ver:   "); Serial.println(sensorGyro.version);
  Serial.print  ("Unique ID:    "); Serial.println(sensorGyro.sensor_id);
  Serial.print  ("Max Value:    "); Serial.print(sensorGyro.max_value); Serial.println(" rad/s");
  Serial.print  ("Min Value:    "); Serial.print(sensorGyro.min_value); Serial.println(" rad/s");
  Serial.print  ("Resolution:   "); Serial.print(sensorGyro.resolution); Serial.println(" rad/s");  
  Serial.println("------------------------------------");
  Serial.println("");
  delay(500);

  //Accel
  sensor_t sensorAccel;
  accel.getSensor(&sensorAccel);
  Serial.println("------------------------------------");
  Serial.print  ("Sensor:       "); Serial.println(sensorAccel.name);
  Serial.print  ("Driver Ver:   "); Serial.println(sensorAccel.version);
  Serial.print  ("Unique ID:    "); Serial.println(sensorAccel.sensor_id);
  Serial.print  ("Max Value:    "); Serial.print(sensorAccel.max_value); Serial.println(" m/s^2");
  Serial.print  ("Min Value:    "); Serial.print(sensorAccel.min_value); Serial.println(" m/s^2");
  Serial.print  ("Resolution:   "); Serial.print(sensorAccel.resolution); Serial.println(" m/s^2");
  Serial.println("------------------------------------");
  Serial.println("");
  delay(500);
}

void setup()
{
  Serial.print("Initializing SD card...");
   pinMode(10, OUTPUT);
 
  if (!SD.begin(10)) {
    Serial.println("initialization failed!");
    return;
  }

  
  //Gyro
    Serial.begin(9600);
    Wire.begin();
    Serial.println("Gyroscope Test");
    
    /* Enable auto-ranging */
    gyro.enableAutoRange(true);
    Serial.println("Test1");
    /* Initialise the sensor */
    if(!gyro.begin())
    {
      /* There was a problem detecting the L3GD20 ... check your connections */
      Serial.println("Ooops, no L3GD20 detected ... Check your wiring!");
      return;
    }
    Serial.println("Test2");
    //Accel
      if(!accel.begin())
      {
        /* There was a problem detecting the ADXL345 ... check your connections */
        Serial.println("Ooops, no LSM303 detected ... Check your wiring!");
        while(1);
      }
    
    //Crap
    /* Display some basic information on this sensor */
    displaySensorDetails();

}



void loop(void)
{
    //Get Data
    sensors_event_t eventGyro; 
    gyro.getEvent(&eventGyro);
    sensors_event_t eventAccel;
    accel.getEvent(&eventAccel);

    char accel1[8];
    char accel2[8];
    char accel3[8];
    char gyro1[8];
    char gyro2[8];
    char gyro3[8];
    dtostrf(eventAccel.acceleration.x,6,3,accel1);
    dtostrf(eventAccel.acceleration.y,6,3,accel2);
    dtostrf(eventAccel.acceleration.z,6,3,accel3);
    dtostrf(eventGyro.gyro.x,6,3,gyro1);
    dtostrf(eventGyro.gyro.y,6,3,gyro2);
    dtostrf(eventGyro.gyro.z,6,3,gyro3);

    snprintf(report, sizeof(report), "Accel(x,y,z:[m/s^2]): %s %s %s   Gyro(x,y,z:[rad/s]): %s %s %s",
    accel1,accel2,accel3,
    gyro1,gyro2,gyro3);
    Serial.println(report);

    
    //Write to SD Card
    myFile = SD.open("IMU.txt", FILE_WRITE);
 
    // if the file opened okay, write to it:
    if (myFile) {
    
      myFile.println(report);
    // close the file:
      myFile.close();
      Serial.println("done.");
    } else {
      // if the file didn't open, print an error:
      Serial.println("error opening file");
    }
    
    delay(500);

}