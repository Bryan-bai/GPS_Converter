// Copyright 2020 DriveX.Tech. All rights reserved.
// 
// Licensed under the License.

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "LatLon_UTM_Converter.hpp"


constexpr char LatLonUTMConverter::ZONE_LETTERS[];


char LatLonUTMConverter::lat_to_zone_letter(double lat)
{
    char zone_letter = -1;

    int lat_d = floor(RAD2DEG(lat));
    if(-80  <= lat_d <= 84) {
        zone_letter = ZONE_LETTERS[(lat_d + 80) / 8];
    }

    return zone_letter;
}

int LatLonUTMConverter::lon_to_zone_number(double lat, double lon)
{
    int zone_number = -1;
    
    int lon_d = floor(RAD2DEG(lon));
    int lat_d = floor(RAD2DEG(lat));

    zone_number = int((lon_d + 180) / 6) + 1;
  
    // special zones
    if(lat_d >= 56 && lat_d < 64 && lon_d >= 3 && lon_d < 12) {
        zone_number = 32;
    } else if(lat_d >= 72 && lat_d <= 84 && lon_d >= 0) {
        if(lon_d < 9) {
            zone_number = 31;
        } else if(lon_d < 21) {
            zone_number = 33;
        } else if(lon_d < 33) {
            zone_number = 35;
        } else if(lon_d < 42) {
            zone_number = 37;
        }
    }
}

static int zone_number_to_central_lon(int zone_number)
{
    assert(zone_number >=1 && zone_number <= 60);
    return (zone_number - 1) * 6 - 180 + 3;  // +3 puts origin in middle of zone
}


void LatLonUTMConverter::LatLon_2_XY(const double lat, const double lon, double& x, double& y, char zone[4])
{
    assert(lat >= DEG2RAD(-80) && lat <= DEG2RAD(84) && lon > -M_PI && lon < M_PI);

    char zone_letter = lat_to_zone_letter(lat);
    int zone_number = lon_to_zone_number(lat, lon);

    double lon_orig = DEG2RAD(zone_number_to_central_lon(zone_number));

    // compute UTM Zone from the latitude and longitude
    sprintf(zone, "%c%d", zone_letter, zone_number);

    double lat_s = sin(lat), lat_c = cos(lat);
    double lat_t = lat_s / lat_c; // tan(lat)
    double lat_t2 = lat_t * lat_t;
    double lat_t4 = lat_t2 * lat_t2;

    double n = R / sqrt(1 - E2 * lat_s * lat_s);
    double c = Ep2 * lat_c * lat_c;

    double a = lat_c * (lon - lon_orig);
    double a2 = a * a;
    double a3 = a * a2;
    double a4 = a * a3;
    double a5 = a * a4;
    double a6 = a * a5;

    double m = R * (M1 * lat - M2 * sin(2.0 * lat) + M3 * sin(4.0 * lat) - M4 * sin(6.0 * lat));
    
    x = 500000.0 + K0 * n *(a + (1.0 - lat_t2 + c) * a3 / 6.0 \
                            + (5.0 - 18.0 * lat_t2 + lat_t4 + 72.0 * c - 58.0 * Ep2) * a5 / 120.0);

    y = K0 * (m + n * lat_t * (a2 / 2.0 + (5.0 - lat_t2 + (9.0 + 4.0 * c) * c) * a4 / 24.0 \
                                + (61.0 - 58.0 * lat_t2 + lat_t4 + 600.0 * c - 330 * Ep2) * a6 / 720.0));

    if(lat < 0) y += 10000000.0; //10000000 meter offset for southern hemisphere
}


void LatLonUTMConverter::XY_2_LatLon(const double x, const double y, const char zone[4], double& lat, double& lon)
{
    char zone_letter;
    int zone_number;

    sscanf(zone, "%c%d", &zone_letter, &zone_number);
    assert(zone_number >=1 && zone_number <= 60);

    double xn = x - 500000.0; //remove 500,000 meter offset from central longitude
    double yn = y;

    if(zone_letter >= 'N') { // in northern hemisphere
    } else { // in southern hemisphere
        yn -= 10000000.0; //remove 10,000,000 meter offset used for southern hemisphere
    }

    double lon_orig = DEG2RAD(zone_number_to_central_lon(zone_number));

    double m = yn / K0;
    double mu = m / (R * M1);

    double p1 = (mu + P2 * sin(2.0 * mu) + P4 * sin(4.0 * mu) + P6 * sin(6.0 * mu) + P8 * sin(8.0 * mu));

    double p1_s = sin(p1), p1_c = cos(p1);
    double p1_s2 = p1_s * p1_s,
            p1_c2 = p1_c * p1_c;
    double p1_t = p1_s / p1_c; // tan(p1)
    double p1_t2 = p1_t * p1_t;
    double p1_t4 = p1_t2 * p1_t2;

    double tmp = 1.0 - E2 * p1_s2;
    double tmp_sqrt = sqrt(tmp);

    double n = R / tmp_sqrt;
    double r = (1.0 - E2) / tmp;
    // double r = R * (1.0 - E2) / (tmp * tmp_sqrt);

    double c = Ep2 * p1_c2;
    double c2 = c * c;

    double d = xn / (n * K0);
    double d2 = d * d;
    double d3 = d * d2;
    double d4 = d * d3;
    double d5 = d * d4;
    double d6 = d * d5;

    lat = p1 - (p1_t / r) * (d2 / 2.0 \
                                - (5.0 + 3.0 * p1_t2 + 10.0 * c - 4.0 * c2 - 9.0 * Ep2) * d4 / 24.0 \
                                + (61.0 + 90.0 * p1_t2 + 298.0 * c + 45.0 * p1_t4 - 252.0 * Ep2 - 3.0 * c2) * d6 / 720.0);
    // lat = p1 - (n * p1_t / r) * (d2 / 2.0 \
    //                                 - (5.0 + 3.0 * p1_t2 + 10.0 * c - 4.0 * c2 - 9.0 * Ep2) * d4 / 24.0 \
    //                                 + (61.0 + 90.0 * p1_t2 + 298.0 * c + 45.0 * p1_t4 - 252.0 * Ep2 - 3.0 * c2) * d6 / 720.0);

    lon = (d - (1.0 + 2.0 * p1_t2 + c) * d3 / 6.0 \
            + (5.0 - 2.0 * c + 28.0 * p1_t2 - 3.0 * c2 + 8.0 * Ep2 + 24.0 * p1_t4) * d5 / 120.0) / p1_c;
    lon += lon_orig;
}
