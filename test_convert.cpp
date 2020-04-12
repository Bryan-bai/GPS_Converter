#include <stdio.h>
#include <iostream>
#include <string>
#include "LatLon_UTM_Converter.hpp"

using namespace std;


int main()
{
    // (5.06132736e+05, 5.81688653e+06)
    double lat = DEG2RAD(52.50207395);
    double lon = DEG2RAD(-2.9096531);

    double x, y;
    char zone[4] = { 0 };
    LatLonUTMConverter::LatLon_2_XY(lat, lon, x, y, zone);

    printf("x: %.8f, y: %.8f, zone: %s\n", x, y, zone);

    double lat_r, lon_r;
    LatLonUTMConverter::XY_2_LatLon(x, y, zone, lat_r, lon_r);
    printf("lat: %.8f, lon: %.8f\n", lat_r, lon_r);
    printf("diff: lat: %.8e, lon: %.8e\n", lat - lat_r, lon - lon_r);

    return 0;
}