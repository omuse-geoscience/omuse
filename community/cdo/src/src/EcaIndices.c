/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

/*
      MODULE      OPERATOR     INDEX    DESCRIPTION
      
      EcaCfd      eca_cfd      CFD      maximum number of consecutive frost days
      EcaCsu      eca_csu      CSU      maximum number of consecutive summer days
      EcaCwdi     eca_cwdi     CWDI     cold wave duration index
      EcaCwfi     eca_cwfi     CWFI     number of cold-spell days
      EcaEtr      eca_etr      ETR      intra-period extreme temperature range
      EcaFd       eca_fd       FD       number of frost days
      EcaGsl      eca_gsl      GSL      growing season length
      EcaHd       eca_hd       HD       heating degree days
      EcaHwdi     eca_hwdi     HWDI     heat wave duration index
      EcaHwfi     eca_hwfi     HWFI     number of warm-spell days
      EcaId       eca_id       ID       number of ice days
      EcaSu       eca_su       SU       number of summer days
      EcaTg10p    eca_tg10p    TG10p    percent of time TX < 10th percentile of daily mean temperature
      EcaTg90p    eca_tg90p    TG90p    percent of time TX > 90th percentile of daily mean temperature
      EcaTn10p    eca_tn10p    TN10p    percent of time TX < 10th percentile of daily minimum temperature   
      EcaTn90p    eca_tn90p    TN90p    percent of time TX > 90th percentile of daily minimum temperature
      EcaTr       eca_tr       TR       number of tropical nights
      EcaTx10p    eca_tx10p    TX10p    percent of time TX < 10th percentile of daily maximum temperature
      EcaTx90p    eca_tx90p    TX90p    percent of time TX > 90th percentile of daily maximum temperature

      EcaCdd      eca_cdd      CDD      maximum number of consecutive dry days
      EcaCwd      eca_cwd      CWD      maximum number of consecutive wet days
      EcaR10mm    eca_r10mm    R10mm    number of days with precipitation >= 10 mm
      EcaR20mm    eca_r20mm    R20mm    number of days with precipitation >= 20 mm
      EcaR75p     eca_r75p     R75p     Percent of time RR > 75th percentile of daily precipitation amount 
      EcaR75ptot  eca_r75ptot  R75pTOT  Percentage of annual total precipitation due to events wit RR > 75th percentile of daily precipitation amount
      EcaR90p     eca_r90p     R90p     Percent of time RR > 90th percentile of daily precipitation amount
      EcaR90ptot  eca_r90ptot  R90pTOT  Percentage of annual total precipitation due to events wit RR > 90th percentile of daily precipitation amount
      EcaR95p     eca_r95p     R95p     Percent of time RR > 95th percentile of daily precipitation amount
      EcaR95ptot  eca_r95ptot  R95pTOT  Percentage of annual total precipitation due to events wit RR > 95th percentile of daily precipitation amount
      EcaR99p     eca_r99p     R99p     Percent of time RR > 75th percentile of daily precipitation amount
      EcaR99ptot  eca_r99ptot  R99pTOT  Percentage of annual total precipitation due to events wit RR > 99th percentile of daily precipitation amount
      EcaRr1      eca_rr1      RR1      number of wet days
      EcaSdii     eca_sdii     SDII     simple daily intensity index
      
      Fdns        fdns                  frost days without surface snow 

      Strwin      strwin                number of strong-wind days
      Strbre      strbre                number of strong-breeze days 
      Strgal      strgal                number of strong-gale days 
      Hurr        hurr                  number of hurricane days 
*/

#include "cdo.h"
#include "cdo_int.h"
#include "ecacore.h"
#include "ecautil.h"


#define TO_DEG_CELSIUS(x) ((x) - 273.15)
#define TO_KELVIN(x) ((x) + 273.15)


static const char CFD_NAME[]         = "consecutive_frost_days_index_per_time_period";
static const char CFD_LONGNAME[]     = "Consecutive frost days index is the greatest number of consecutive frost days in a given time period. Frost days is the number of days where minimum of temperature is below 0 degree Celsius. The time period should be defined by the bounds of the time coordinate.";
static const char CFD_UNITS[]        = "No.";
static const char CFD_NAME2[]        = "number_of_cfd_periods_with_more_than_5days_per_time_period";
static const char CFD_LONGNAME2[]    = "Number of cfd periods in given time period with more than 5 days. The time period should be defined by the bounds of the time coordinate.";
static const char CFD_UNITS2[]       = "No.";

static const char CSU_NAME[]         = "consecutive_summer_days_index_per_time_period";
static const char CSU_LONGNAME[]     = "Consecutive summer days index is the greatest number of consecutive summer days in a given time period. Summer days is the number of days where maximum of temperature is above 25 degree Celsius. The time period should be defined by the bounds of the time coordinate.";
static const char CSU_UNITS[]        = "No.";
static const char CSU_NAME2[]        = "number_of_csu_periods_with_more_than_5days_per_time_period";
static const char CSU_LONGNAME2[]    = "Number of csu periods in given time period with more than 5 days. The time period should be defined by the bounds of the time coordinate.";
static const char CSU_UNITS2[]       = "No.";

static const char CWDI_NAME[]        = "cold_wave_duration_index_wrt_mean_of_reference_period";
static const char CWDI_LONGNAME[]    = "This is the number of days per time period where in intervals of at least %d consecutive days the daily minimum temperature is more than %1.0f degrees below a reference value. The reference value is calculated  as the mean of minimum temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char CWDI_UNITS[]       = "No.";
static const char CWDI_NAME2[]       = "cold_waves_per_time_period";
static const char CWDI_LONGNAME2[]   = "Number of cold waves per time period. The time period should be defined by the bounds of the time coordinate.";
static const char CWDI_UNITS2[]      = "No.";

static const char CWFI_NAME[]        = "cold_spell_days_index_wrt_10th_percentile_of_reference_period";
static const char CWFI_LONGNAME[]    = "This is the number of days per time period where in intervals of at least %d consecutive days the daily mean temperature is below a reference value. The reference value is calculated  as the 10th percentile of daily mean temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char CWFI_UNITS[]       = "No.";
static const char CWFI_NAME2[]       = "cold_spell_periods_per_time_period";
static const char CWFI_LONGNAME2[]   = "Number of cold spell periods per time period. The time period should be defined by the bounds of the time coordinate.";
static const char CWFI_UNITS2[]      = "No.";

static const char ETR_NAME[]         = "intra_period_extreme_temperature_range";
static const char ETR_LONGNAME[]     = "Difference between the absolute extreme temperatures in observation period. The time period should be defined by the bounds of the time coordinate.";
static const char ETR_UNITS[]        = "K";

static const char FD_NAME[]          = "frost_days_index_per_time_period";
static const char FD_LONGNAME[]      = "Frost days index is the number of days where minimum of temperature is below 0 degree Celsius. The time period should be defined by the bounds of the time coordinate.";
static const char FD_UNITS[]         = "No.";

static const char GSL_NAME[]         = "thermal_growing_season_length";
static const char GSL_LONGNAME[]     = "Counted are the number of days per calendar year between the first occurrence of at least %d consecutive days where the daily mean temperature is above %1.0f degree Celsius and the first occurrence of at least %d consecutive days after 1st of July where the daily mean temperature is below %1.0f degree Celsius. The time period should be defined by the bounds of the time coordinate.";
static const char GSL_UNITS[]        = "No.";
static const char GSL_NAME2[]        = "day_of_year_of_growing_season_start";
static const char GSL_LONGNAME2[]    = "Day of year of growing season start. The time period should be defined by the bounds of the time coordinate.";
static const char GSL_UNITS2[]       = "No.";

static const char HD_NAME[]          = "heating_degree_days_per_time_period";
static const char HD_LONGNAME[]      = "Heating degree days relates the outside temperature with the room temperature during the heating period. It is the sum of the difference between room temperature X and daily mean temperature Y on days where Y is below a given constant A. X is 20 degree Celsius and A is 15 degree Celsius according to VDI guidelines. According to ECAD both X and A are 17 degree Celsius. The time period should be defined by the bounds of the time coordinate.";
static const char HD_UNITS[]         = "No.";

static const char HWDI_NAME[]        = "heat_wave_duration_index_wrt_mean_of_reference_period";
static const char HWDI_LONGNAME[]    = "This is the number of days per time period where in intervals of at least %d consecutive days the daily maximum temperature is more than %1.0f degrees above a reference value. The reference value is calculated  as the mean of maximum temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char HWDI_UNITS[]       = "No.";
static const char HWDI_NAME2[]       = "heat_waves_per_time_period";
static const char HWDI_LONGNAME2[]   = "Number of heat waves per time period. The time period should be defined by the bounds of the time coordinate.";
static const char HWDI_UNITS2[]      = "No.";

static const char HWFI_NAME[]        = "warm_spell_days_index_wrt_90th_percentile_of_reference_period";
static const char HWFI_LONGNAME[]    = "This is the number of days per time period where in intervals of at least %d consecutive days the daily mean temperature is above a reference value. The reference value is calculated  as the 90th percentile of daily mean temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char HWFI_UNITS[]       = "No.";
static const char HWFI_NAME2[]       = "warm_spell_periods_per_time_period";
static const char HWFI_LONGNAME2[]   = "Number of warm spell periods per time period. The time period should be defined by the bounds of the time coordinate.";
static const char HWFI_UNITS2[]      = "No.";

static const char ID_NAME[]          = "ice_days_index_per_time_period";
static const char ID_LONGNAME[]      = "Ice days index is the number of days where maximum of temperature is below 0 degree Celsius. The time period should be defined by the bounds of the time coordinate.";
static const char ID_UNITS[]         = "No.";

static const char SU_NAME[]          = "summer_days_index_per_time_period";
static const char SU_LONGNAME[]      = "Summer days index is the number of days where maximum of temperature is above %1.0f degree Celsius. The time period should be defined by the bounds of the time coordinate.";
static const char SU_UNITS[]         = "No.";

static const char TG10P_NAME[]       = "cold_days_percent_wrt_10th_percentile_of_reference_period";
static const char TG10P_LONGNAME[]   = "This is the percent of time per time period where daily mean temperature is below a reference value. The reference value is calculated  as the 10th percentile of daily mean temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char TG10P_UNITS[]      = "Percent";

static const char TG90P_NAME[]       = "warm_days_percent_wrt_90th_percentile_of_reference_period";
static const char TG90P_LONGNAME[]   = "This is the percent of time per time period where daily mean  temperature is above a reference value. The reference value is calculated  as the 90th percentile of daily mean temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char TG90P_UNITS[]      = "Percent";

static const char TN10P_NAME[]       = "cold_nights_percent_wrt_10th_percentile_of_reference_period";
static const char TN10P_LONGNAME[]   = "This is the percent of time per time period where daily minimum  temperature is below a reference value. The reference value is calculated  as the 10th percentile of daily minimum temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char TN10P_UNITS[]      = "Percent";

static const char TN90P_NAME[]       = "warm_nights_percent_wrt_90th_percentile_of_reference_period";
static const char TN90P_LONGNAME[]   = "This is the percent of time per time period where daily minimum  temperature is above a reference value. The reference value is calculated  as the 90th percentile of daily minimum temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char TN90P_UNITS[]      = "Percent";

static const char TR_NAME[]          = "tropical_nights_index_per_time_period";
static const char TR_LONGNAME[]      = "Tropical nights index is the number of days where minimum of temperature is above %1.0f degree Celsius. The time period should be defined by the bounds of the time coordinate.";
static const char TR_UNITS[]         = "No.";

static const char TX10P_NAME[]       = "very_cold_days_percent_wrt_10th_percentile_of_reference_period";
static const char TX10P_LONGNAME[]   = "This is the percent of time per time period where daily maximum temperature is below a reference value. The reference value is calculated  as the 10th percentile of daily maximum temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char TX10P_UNITS[]      = "Percent";

static const char TX90P_NAME[]       = "very_warm_days_percent_wrt_90th_percentile_of_reference_period";
static const char TX90P_LONGNAME[]   = "This is the percent of time per time period where daily maximum  temperature is above a reference value. The reference value is calculated  as the 90th percentile of daily maximum temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char TX90P_UNITS[]      = "Percent";

static const char CDD_NAME[]         = "consecutive_dry_days_index_per_time_period";
static const char CDD_LONGNAME[]     = "Consecutive dry days is the greatest number of consecutive days per time period with daily precipitation amount  below %g mm. The time period should be defined by the bounds of the time coordinate.";
static const char CDD_UNITS[]        = "No.";
static const char CDD_NAME2[]        = "number_of_cdd_periods_with_more_than_5days_per_time_period";
static const char CDD_LONGNAME2[]    = "Number of cdd periods in given time period with more than 5 days. The time period should be defined by the bounds of the time coordinate.";
static const char CDD_UNITS2[]       = "No.";

static const char CWD_NAME[]         = "consecutive_wet_days_index_per_time_period";
static const char CWD_LONGNAME[]     = "Consecutive wet days is the greatest number of consecutive days per time period with daily precipitation above %g mm. The time period should be defined by the bounds of the time coordinate.";
static const char CWD_UNITS[]        = "No.";
static const char CWD_NAME2[]        = "number_of_cwd_periods_with_more_than_5days_per_time_period";
static const char CWD_LONGNAME2[]    = "Number of cwd periods in given time period with more than 5 days. The time period should be defined by the bounds of the time coordinate.";
static const char CWD_UNITS2[]       = "No.";

static const char PD_NAME[]          = "precipitation_days_index_per_time_period";
static const char PD_LONGNAME[]      = "precipitation days is the number of days per time period with daily precipitation sum exceeding %g mm. The time period should be defined by the bounds of the time coordinate.";
static const char PD_UNITS[]         = "No.";

static const char R10MM_NAME[]       = "heavy_precipitation_days_index_per_time_period";
static const char R10MM_LONGNAME[]   = "Heavy precipitation days is the number of days per time period with daily precipitation sum exceeding 10 mm. The time period should be defined by the bounds of the time coordinate.";
static const char R10MM_UNITS[]      = "No.";

static const char R20MM_NAME[]       = "very_heavy_precipitation_days_index_per_time_period";
static const char R20MM_LONGNAME[]   = "Very heavy precipitation days is the number of days with daily precipitation sum exceeding 20 mm. The time period should be defined by the bounds of the time coordinate.";
static const char R20MM_UNITS[]      = "No.";

static const char R75P_NAME[]        = "moderate_wet_days_wrt_75th_percentile_of_reference_period";
static const char R75P_LONGNAME[]    = "This is the percent of time per time period of wet days (daily sum at least 1 mm / day) where daily precipitation amount of a wet day is above a reference value. The reference value is calculated  as the 75th percentile of all wet days of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char R75P_UNITS[]       = "Percent";

static const char R75PTOT_NAME[]     = "precipitation_percent_due_to_R75p_days";
static const char R75PTOT_LONGNAME[] = "Percentage of  total precipitation amount per time period due to moderate_wet_days_wrt_75th_percentile_of_reference_period. The time period should be defined by the bounds of the time coordinate.";
static const char R75PTOT_UNITS[]    = "Percent";

static const char R90P_NAME[]        = "wet_days_wrt_90th_percentile_of_reference_period";
static const char R90P_LONGNAME[]    = "This is the percent of time per time period of wet days (daily sum at least 1 mm / day) where daily precipitation amount of a wet day is above a reference value. The reference value is calculated  as the 90th percentile of all wet days of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char R90P_UNITS[]       = "Percent";

static const char R90PTOT_NAME[]     = "precipitation_percent_due_to_R90p_days";
static const char R90PTOT_LONGNAME[] = "Percentage of  total precipitation  amount per time period  due towet_days_wrt_90th_percentile_of_reference_period. The time period should be defined by the bounds of the time coordinate.";
static const char R90PTOT_UNITS[]    = "Percent";

static const char R95P_NAME[]        = "very_wet_days_wrt_95th_percentile_of_reference_period";
static const char R95P_LONGNAME[]    = "This is the percent of time per time period of wet days (daily sum at least 1 mm / day) where daily precipitation amount of a wet day is above a reference value. The reference value is calculated  as the 95th percentile of all wet days of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char R95P_UNITS[]       = "Percent";

static const char R95PTOT_NAME[]     = "precipitation_percent_due_to_R95p_days";
static const char R95PTOT_LONGNAME[] = "Percentage of  total  precipitation amount per time period  due to  very_wet_days_wrt_95th_percentile_of_reference_period. The time period should be defined by the bounds of the time coordinate.";
static const char R95PTOT_UNITS[]    = "Percent";

static const char R99P_NAME[]        = "extremely_wet_days_wrt_99th_percentile_of_reference_period";
static const char R99P_LONGNAME[]    = "This is the percent of time per time period of wet days (daily sum at least 1 mm / day) where daily precipitation amount of a wet day is above a reference value. The reference value is calculated  as the 99th percentile of all wet days of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char R99P_UNITS[]       = "Percent";

static const char R99PTOT_NAME[]     = "precipitation_percent_due_to_R99p_days";
static const char R99PTOT_LONGNAME[] = "percentage of  total  precipitation amount per time period  due to  extremely_wet_days_wrt_99th_percentile_of_reference_period. The time period should be defined by the bounds of the time coordinate.";
static const char R99PTOT_UNITS[]    = "Percent";

static const char RR1_NAME[]         = "wet_days_index_per_time_period";
static const char RR1_LONGNAME[]     = "Wet days index is the number of days per time period with daily precipitation of at least %g mm.  The time period should be defined by the bounds of the time coordinate.";
static const char RR1_UNITS[]        = "No.";

static const char RX1DAY_NAME[]      = "highest_one_day_precipitation_amount_per_time_period";
static const char RX1DAY_LONGNAME[]  = "Highest one  day precipitation  is the maximum of one day precipitation amount in a given time period. The time period should be defined by the bounds of the time coordinate.";
static const char RX1DAY_UNITS[]     = "mm per day";

static const char RX5DAY_NAME[]      = "highest_five_day_precipitation_amount_per_time_period";
static const char RX5DAY_LONGNAME[]  = "Highest precipitation amount  for five day interval (including the calendar day as the last day). The time period should be defined by the bounds of the time coordinate.";
static const char RX5DAY_UNITS[]     = "mm per 5 day";
static const char RX5DAY_NAME2[]     = "number_of_5day_heavy_precipitation_periods_per_time_period";
static const char RX5DAY_LONGNAME2[] = "Number of 5day periods in given time period with precipitation amount exceeding %1.0f mm / 5 days. The time period should be defined by the bounds of the time coordinate.";
static const char RX5DAY_UNITS2[]    = "No.";

static const char SDII_NAME[]        = "simple_daily_intensitiy_index_per_time_period";
static const char SDII_LONGNAME[]    = "Simple daily intensity index is the mean of precipitation amount on wet days. A wet day is a day with precipitation sum of at least %g mm. The time period should be defined by the bounds of the time coordinate.";
static const char SDII_UNITS[]       = "mm";

static const char FDNS_NAME[]        = "frost_days_where_no_snow_index_per_time_period";       
static const char FDNS_LONGNAME[]    = "Frost days where no snow index is the number of days without snowcover and where the minimum of temperature is below 0 degree Celsius. The time period should be defined by the bounds of the time coordinate.";       
static const char FDNS_UNITS[]       = "No.";

static const char STRWIN_NAME[]      = "strong_wind_days_index_per_time_period";
static const char STRWIN_LONGNAME[]  = "Strong wind days index is the number of days per time period where maximum wind speed is above %1.0f m/s. The time period should be defined by the bounds of the time coordinate.";
static const char STRWIN_UNITS[]     = "No.";
static const char STRWIN_NAME2[]     = "consecutive_strong_wind_days_index_per_time_period";
static const char STRWIN_LONGNAME2[] = "Greatest number of consecutive strong wind days per time period. The time period should be defined by the bounds of the time coordinate.";
static const char STRWIN_UNITS2[]    = "No.";

static const char STRBRE_NAME[]      = "strong_breeze_days_index_per_time_period";
static const char STRBRE_LONGNAME[]  = "Strong breeze days index is the number of days per time period where maximum wind speed is above 10.5 m/s. The time period should be defined by the bounds of the time coordinate.";
static const char STRBRE_NAME2[]     = "consecutive_strong_breeze_days_index_per_time_period";
static const char STRBRE_LONGNAME2[] = "Greatest number of consecutive strong breeze days per time period. The time period should be defined by the bounds of the time coordinate.";

static const char STRGAL_NAME[]      = "strong_gale_days_index_per_time_period";
static const char STRGAL_LONGNAME[]  = "Strong gale days index is the number of days per time period where maximum wind speed is above 20.5 m/s. The time period should be defined by the bounds of the time coordinate.";
static const char STRGAL_NAME2[]     = "consecutive_strong_gale_days_index_per_time_period";
static const char STRGAL_LONGNAME2[] = "Greatest number of consecutive strong gale days per time period. The time period should be defined by the bounds of the time coordinate.";

static const char HURR_NAME[]        = "hurricane_days_index_per_time_period";
static const char HURR_LONGNAME[]    = "Hurricane days index is the number of days per time period where maximum wind speed is above 32.5 m/s. The time period should be defined by the bounds of the time coordinate.";
static const char HURR_NAME2[]       = "consecutive_hurricane_days_index_per_time_period";
static const char HURR_LONGNAME2[]   = "Greatest number of consecutive hurricane days per time period. The time period should be defined by the bounds of the time coordinate.";


/* ECA temperature indices */


void *EcaCfd(void *argument)
{
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_cfd", 0, 31, NULL);
  
  request.var1.name     = CFD_NAME;
  request.var1.longname = CFD_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = farselltc;
  request.var1.f1arg    = TO_KELVIN(0.0);
  request.var1.f2       = farnum2;
  request.var1.f3       = farmax;
  request.var1.mulc     = 0.0;
  request.var1.addc     = 0.0;
  request.var1.epilog   = NONE;
  request.var2.name     = CFD_NAME2;
  request.var2.longname = CFD_LONGNAME2;
  request.var2.units    = CFD_UNITS2;
  request.var2.h1       = farseleqc;
  request.var2.h1arg    = 6;
  request.var2.h2       = NULL;
  request.var2.h3       = farnum;
   
  eca1(&request);
  cdoFinish();
  
  return (0);
}


void *EcaCsu(void *argument)
{
  double argT = 25.0;
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_csu", 0, 31, NULL);

  if ( operatorArgc() > 0 ) argT = parameter2double(operatorArgv()[0]);

  request.var1.name     = CSU_NAME;
  request.var1.longname = CSU_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = farselgtc;
  request.var1.f1arg    = TO_KELVIN(argT);
  request.var1.f2       = farnum2;
  request.var1.f3       = farmax;
  request.var1.mulc     = 0.0;
  request.var1.addc     = 0.0;
  request.var1.epilog   = NONE;
  request.var2.name     = CSU_NAME2;
  request.var2.longname = CSU_LONGNAME2;
  request.var2.units    = CSU_UNITS2;
  request.var2.h1       = farseleqc;
  request.var2.h1arg    = 6;
  request.var2.h2       = NULL;
  request.var2.h3       = farnum;
  
  eca1(&request);
  cdoFinish();
  
  return (0);
}


void *EcaCwdi(void *argument)
{
  char *longname;
  int argN = 6;
  double argT = 5.0;
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_cwdi", 0, 31, NULL);

  if ( operatorArgc() > 0 ) argN = parameter2int(operatorArgv()[0]);
  if ( operatorArgc() > 1 ) argT = parameter2double(operatorArgv()[1]);
  
  longname = (char*) malloc(strlen(CWDI_LONGNAME) + 80);
  sprintf(longname, CWDI_LONGNAME, argN, argT);

  request.var1.name     = CWDI_NAME;
  request.var1.longname = longname;
  request.var1.units    = CWDI_UNITS;
  request.var1.f1       = NULL;
  request.var1.f2       = farcsub;
  request.var1.f2arg    = argT;
  request.var1.f3       = farsellt;
  request.var1.f4       = farnum2;
  request.var1.f5       = farnum3;
  request.var1.f5arg    = argN;
  request.var1.epilog   = NONE;
  request.var2.name     = CWDI_NAME2;
  request.var2.longname = CWDI_LONGNAME2;
  request.var2.units    = CWDI_UNITS2;
  request.var2.h1       = farseleqc;
  request.var2.h1arg    = argN;
  request.var2.h2       = farnum;
   
  eca2(&request);
  
  free(longname);
  cdoFinish();
  
  return (0);
}


void *EcaCwfi(void *argument)
{
  char *longname;
  int argN = 6;
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_cwfi", 0, 31, NULL);

  if ( operatorArgc() > 0 ) argN = parameter2int(operatorArgv()[0]);

  longname = (char*) malloc(strlen(CWFI_LONGNAME) + 40);
  sprintf(longname, CWFI_LONGNAME, argN);

  request.var1.name     = CWFI_NAME;
  request.var1.longname = longname;
  request.var1.units    = CWFI_UNITS;
  request.var1.f1       = NULL;
  request.var1.f2       = NULL;
  request.var1.f3       = farsellt;
  request.var1.f4       = farnum2;
  request.var1.f5       = farnum3;
  request.var1.f5arg    = argN;
  request.var1.epilog   = NONE;
  request.var2.name     = CWFI_NAME2;
  request.var2.longname = CWFI_LONGNAME2;
  request.var2.units    = CWFI_UNITS2;
  request.var2.h1       = farseleqc;
  request.var2.h1arg    = argN;
  request.var2.h2       = farnum;
   
  eca2(&request);
  
  free(longname);
  cdoFinish();
  
  return (0);
}


void *EcaEtr(void *argument)
{
  ECA_REQUEST_3 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_etr", 0, 31, NULL);
  
  request.name     = ETR_NAME;
  request.longname = ETR_LONGNAME;
  request.units    = NULL;
  request.f1       = farmax; 
  request.f2       = farmin;
  request.f3       = farsub;
   
  eca3(&request);
  cdoFinish();
  
  return (0);
}


void *EcaFd(void *argument)
{
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_fd", 0, 31, NULL);
  
  request.var1.name     = FD_NAME;
  request.var1.longname = FD_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = farselltc; 
  request.var1.f1arg    = TO_KELVIN(0.0);
  request.var1.f2       = farnum;
  request.var1.f3       = NULL;
  request.var1.mulc     = 0.0;    
  request.var1.addc     = 0.0;
  request.var1.epilog   = NONE;    
  request.var2.h2       = NULL; 
  request.var2.h3       = NULL; 
   
  eca1(&request);
  cdoFinish();
  
  return (0);
}


/*
 * Definition of GSL: (Thermal) Growing Season Length start at the first span
 * of at least 6 (argN) days with T > 5.0°C (argT) in first half of the year
 * and ends at the first span of ar least 6 (argN) days with T < 5.0°C (argT)
 * in the second half.
 * ATTENTION: Year of the northern hemisphere starts in january to
 * december, whereas for the southern hemisphere is goes from july to june!
 * Hence, at least 18 Month of data is needed for computing the gsl of the
 * whole earth.
*/
void *EcaGsl(void *argument)
{
  char *longname;
  int argN = 6; 
  double argT = 5.0;
  double minLandFraction = 0.5;
  ECA_REQUEST_4 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_gsl", 0, 10, NULL);
  
  if ( operatorArgc() > 0 ) argN = parameter2int(operatorArgv()[0]);
  if ( operatorArgc() > 1 ) argT = parameter2double(operatorArgv()[1]);
  if ( operatorArgc() > 2 ) minLandFraction = parameter2double(operatorArgv()[2]);

  longname = (char*) malloc(strlen(GSL_LONGNAME) + 160);
  sprintf(longname, GSL_LONGNAME, argN, argT, argN, argT);
  
  request.name      = GSL_NAME;
  request.longname  = longname;
  request.units     = GSL_UNITS;
  request.name2     = GSL_NAME2;
  request.longname2 = GSL_LONGNAME2;
  request.units2    = GSL_UNITS2;
  request.s1        = farselgtc; 
  request.s1arg     = TO_KELVIN(argT);
  request.s2        = farselltc;
  request.s2arg     = TO_KELVIN(argT);
  request.s3        = farselgec;
  request.s3arg     = minLandFraction;
  request.consecutiveDays = argN;    
   
  eca4(&request);
  
  free(longname);

  cdoFinish();
  
  return (0);
}


void *EcaHd(void *argument)
{
  double argX = 17.0;
  double argA = 17.0;
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_hd", 0, 31, NULL);

  if ( operatorArgc() > 0 ) 
    {
      argX = parameter2double(operatorArgv()[0]);
      argA = argX;
    }
  if ( operatorArgc() > 1 ) 
    argA = parameter2double(operatorArgv()[1]);
  
  request.var1.name     = HD_NAME;
  request.var1.longname = HD_LONGNAME;
  request.var1.units    = HD_UNITS;
  request.var1.f1       = farselltc; 
  request.var1.f1arg    = TO_KELVIN(argA);
  request.var1.f2       = farsum;
  request.var1.f3       = NULL;
  request.var1.mulc     = -1.0;    
  request.var1.addc     = TO_KELVIN(argX);
  request.var1.epilog   = NONE;    
  request.var2.h2       = NULL; 
  request.var2.h3       = NULL; 
   
  eca1(&request);

  cdoFinish();
  
  return (0);
}


void *EcaHwdi(void *argument)
{
  char *longname;
  int argN = 6;
  double argT = 5.0;
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_hwdi", 0, 31, NULL);

  if ( operatorArgc() > 0 ) argN = parameter2int(operatorArgv()[0]);
  if ( operatorArgc() > 1 ) argT = parameter2double(operatorArgv()[1]);
  
  longname = (char*) malloc(strlen(HWDI_LONGNAME) + 80);
  sprintf(longname, HWDI_LONGNAME, argN, argT);
  
  request.var1.name     = HWDI_NAME;
  request.var1.longname = longname;
  request.var1.units    = HWDI_UNITS;
  request.var1.f1       = NULL;
  request.var1.f2       = farcadd;
  request.var1.f2arg    = argT;
  request.var1.f3       = farselgt;
  request.var1.f4       = farnum2;
  request.var1.f5       = farnum3;
  request.var1.f5arg    = argN;
  request.var1.epilog   = NONE;
  request.var2.name     = HWDI_NAME2;
  request.var2.longname = HWDI_LONGNAME2;
  request.var2.units    = HWDI_UNITS2;
  request.var2.h1       = farseleqc;
  request.var2.h1arg    = argN;
  request.var2.h2       = farnum;
   
  eca2(&request);
  
  free(longname);
  cdoFinish();
  
  return (0);
}


void *EcaHwfi(void *argument)
{
  char *longname;
  int argN = 6;
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_hwfi", 0, 31, NULL);

  if ( operatorArgc() > 0 ) argN = parameter2int(operatorArgv()[0]);

  longname = (char*) malloc(strlen(HWFI_LONGNAME) + 40);
  sprintf(longname, HWFI_LONGNAME, argN);

  request.var1.name     = HWFI_NAME;
  request.var1.longname = longname;
  request.var1.units    = HWFI_UNITS;
  request.var1.f1       = NULL;
  request.var1.f2       = NULL;
  request.var1.f3       = farselgt;
  request.var1.f4       = farnum2;
  request.var1.f5       = farnum3;
  request.var1.f5arg    = argN;
  request.var1.epilog   = NONE;
  request.var2.name     = HWFI_NAME2;
  request.var2.longname = HWFI_LONGNAME2;
  request.var2.units    = HWFI_UNITS2;
  request.var2.h1       = farseleqc;
  request.var2.h1arg    = argN;
  request.var2.h2       = farnum;
   
  eca2(&request);
  
  free(longname);
  cdoFinish();
  
  return (0);
}


void *EcaId(void *argument)
{
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_id", 0, 31, NULL);

  request.var1.name     = ID_NAME;
  request.var1.longname = ID_LONGNAME;
  request.var1.units    = ID_UNITS;
  request.var1.f1       = farselltc;
  request.var1.f1arg    = TO_KELVIN(0.0);
  request.var1.f2       = farnum;
  request.var1.f3       = NULL;
  request.var1.mulc     = 0.0;    
  request.var1.addc     = 0.0;   
  request.var1.epilog   = NONE; 
  request.var2.h2       = NULL; 
  request.var2.h3       = NULL; 
    
  eca1(&request);
  cdoFinish();
  
  return (0);
}


void *EcaSu(void *argument)
{
  char *longname;
  double argT = 25.0;
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_su", 0, 31, NULL);

  if ( operatorArgc() > 0 ) argT = parameter2double(operatorArgv()[0]);
  longname = (char*) malloc(strlen(SU_LONGNAME) + 40);
  sprintf(longname, SU_LONGNAME, argT);

  request.var1.name     = SU_NAME;
  request.var1.longname = longname;
  request.var1.units    = NULL;
  request.var1.f1       = farselgtc;
  request.var1.f1arg    = TO_KELVIN(argT);
  request.var1.f2       = farnum;
  request.var1.f3       = NULL;
  request.var1.mulc     = 0.0;    
  request.var1.addc     = 0.0;    
  request.var1.epilog   = NONE; 
  request.var2.h2       = NULL; 
  request.var2.h3       = NULL; 
 
  eca1(&request);
  
  free(longname);
  cdoFinish();
  
  return (0);
}


void *EcaTg10p(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_tg10p", 0, 31, NULL);

  request.var1.name     = TG10P_NAME;
  request.var1.longname = TG10P_LONGNAME;
  request.var1.units    = TG10P_UNITS;
  request.var1.f1       = NULL;
  request.var1.f2       = NULL;
  request.var1.f3       = farsellt;
  request.var1.f4       = farnum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TIME;
  request.var2.h2       = NULL;
    
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaTg90p(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_tg90p", 0, 31, NULL);

  request.var1.name     = TG90P_NAME;
  request.var1.longname = TG90P_LONGNAME;
  request.var1.units    = TG90P_UNITS;
  request.var1.f1       = NULL;
  request.var1.f2       = NULL;
  request.var1.f3       = farselgt;
  request.var1.f4       = farnum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TIME;
  request.var2.h2       = NULL;
  
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaTn10p(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_tn10p", 0, 31, NULL);

  request.var1.name     = TN10P_NAME;
  request.var1.longname = TN10P_LONGNAME;
  request.var1.units    = TN10P_UNITS;
  request.var1.f1       = NULL;
  request.var1.f2       = NULL;
  request.var1.f3       = farsellt;
  request.var1.f4       = farnum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TIME;
  request.var2.h2       = NULL;
    
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaTn90p(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_tn90p", 0, 31, NULL);

  request.var1.name     = TN90P_NAME;
  request.var1.longname = TN90P_LONGNAME;
  request.var1.units    = TN90P_UNITS;
  request.var1.f1       = NULL;
  request.var1.f2       = NULL;
  request.var1.f3       = farselgt;
  request.var1.f4       = farnum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TIME;
  request.var2.h2       = NULL;
    
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaTr(void *argument)
{
  char *longname;
  double argT = 20.0;
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_tr", 0, 31, NULL);

  if ( operatorArgc() > 0 ) argT = parameter2double(operatorArgv()[0]);
  longname = (char*) malloc(strlen(TR_LONGNAME) + 40);
  sprintf(longname, TR_LONGNAME, argT);
 
  request.var1.name     = TR_NAME;
  request.var1.longname = longname;
  request.var1.units    = TR_UNITS;
  request.var1.f1       = farselgtc;
  request.var1.f1arg    = TO_KELVIN(argT);
  request.var1.f2       = farnum;
  request.var1.f3       = NULL;
  request.var1.mulc     = 0.0;    
  request.var1.addc     = 0.0;    
  request.var1.epilog   = NONE; 
  request.var2.h2       = NULL; 
  request.var2.h3       = NULL; 
   
  eca1(&request);
  
  free(longname);
  cdoFinish();
  
  return (0);
}


void *EcaTx10p(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_tx10p", 0, 31, NULL);

  request.var1.name     = TX10P_NAME;
  request.var1.longname = TX10P_LONGNAME;
  request.var1.units    = TX10P_UNITS;
  request.var1.f1       = NULL;
  request.var1.f2       = NULL;
  request.var1.f3       = farsellt;
  request.var1.f4       = farnum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TIME;
  request.var2.h2       = NULL;
    
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaTx90p(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_tx90p", 0, 31, NULL);
 
  request.var1.name     = TX90P_NAME;
  request.var1.longname = TX90P_LONGNAME;
  request.var1.units    = TX90P_UNITS;
  request.var1.f1       = NULL;
  request.var1.f2       = NULL;
  request.var1.f3       = farselgt;
  request.var1.f4       = farnum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TIME;
  request.var2.h2       = NULL;
   
  eca2(&request);
  cdoFinish();
  
  return (0);
}


/* ECA precipitation indices */


void *EcaCdd(void *argument)
{
  ECA_REQUEST_1 request;
  char lnamebuffer[1024];
  double threshold = 1;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_cdd", 0, 31, NULL);

  if ( operatorArgc() == 1 ) threshold = parameter2double(operatorArgv()[0]);
  if ( operatorArgc() > 1 ) cdoAbort("Too many arguments!");

  sprintf(lnamebuffer, CDD_LONGNAME, threshold);

  request.var1.name     = CDD_NAME;
  request.var1.longname = lnamebuffer;
  request.var1.units    = CDD_UNITS;
  request.var1.f1       = farselltc;
  request.var1.f1arg    = threshold;
  request.var1.f2       = farnum2;
  request.var1.f3       = farmax;
  request.var1.mulc     = 0.0;
  request.var1.addc     = 0.0;
  request.var1.epilog   = NONE;
  request.var2.name     = CDD_NAME2;
  request.var2.longname = CDD_LONGNAME2;
  request.var2.units    = CDD_UNITS2;
  request.var2.h1       = farseleqc;
  request.var2.h1arg    = 6;
  request.var2.h2       = NULL;
  request.var2.h3       = farnum;
   
  eca1(&request);
  cdoFinish();
  
  return (0);
}


void *EcaCwd(void *argument)
{
  ECA_REQUEST_1 request;
  char lnamebuffer[1024];
  double threshold = 1;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_cwd", 0, 31, NULL);

  if ( operatorArgc() == 1 ) threshold = parameter2double(operatorArgv()[0]);
  if ( operatorArgc() > 1 ) cdoAbort("Too many arguments!");

  sprintf(lnamebuffer, CWD_LONGNAME, threshold);

  request.var1.name     = CWD_NAME;
  request.var1.longname = lnamebuffer;
  request.var1.units    = CWD_UNITS;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = threshold;
  request.var1.f2       = farnum2;
  request.var1.f3       = farmax;
  request.var1.mulc     = 0.0;
  request.var1.addc     = 0.0;
  request.var1.epilog   = NONE;
  request.var2.name     = CWD_NAME2;
  request.var2.longname = CWD_LONGNAME2;
  request.var2.units    = CWD_UNITS2;
  request.var2.h1       = farseleqc;
  request.var2.h1arg    = 6;
  request.var2.h2       = NULL;
  request.var2.h3       = farnum;
   
  eca1(&request);
  cdoFinish();
  
  return (0);
}


void *EcaPd(void *argument)
{
  int ECA_PD, ECA_R10MM, ECA_R20MM;
  int operatorID;
  char lnamebuffer[1024];
  double threshold = 0;
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);

  ECA_PD      = cdoOperatorAdd("eca_pd",      0, 31, NULL);
  ECA_R10MM   = cdoOperatorAdd("eca_r10mm",   0, 31, NULL);
  ECA_R20MM   = cdoOperatorAdd("eca_r20mm",   0, 31, NULL);

  operatorID = cdoOperatorID();

  if ( operatorID == ECA_PD )
    {
      operatorInputArg("daily precipitation amount threshold in [mm]");

      if ( operatorArgc() < 1 ) cdoAbort("Too few arguments!");
      if ( operatorArgc() > 1 ) cdoAbort("Too many arguments!");

      threshold = parameter2double(operatorArgv()[0]);

      if ( threshold < 0 ) cdoAbort("Parameter out of range: threshold = %d", threshold);

      sprintf(lnamebuffer, PD_LONGNAME, threshold);
      request.var1.name     = PD_NAME;
      request.var1.longname = lnamebuffer;
      request.var1.units    = PD_UNITS;
    }
  else if ( operatorID == ECA_R10MM )
    {
      threshold = 10;

      request.var1.name     = R10MM_NAME;
      request.var1.longname = R10MM_LONGNAME;
      request.var1.units    = R10MM_UNITS;
    }
  else if ( operatorID == ECA_R20MM )
    {
      threshold = 20;

      request.var1.name     = R20MM_NAME;
      request.var1.longname = R20MM_LONGNAME;
      request.var1.units    = R20MM_UNITS;
    }
  
  if ( cdoVerbose ) cdoPrint("threshold = %g", threshold);

  request.var1.f1       = farselgec;
  request.var1.f1arg    = threshold;
  request.var1.f2       = farnum;
  request.var1.f3       = NULL;
  request.var1.mulc     = 0.0;   
  request.var1.addc     = 0.0;    
  request.var1.epilog   = NONE;
  request.var2.h2       = NULL;
  request.var2.h3       = NULL;
  
  eca1(&request);
  cdoFinish();
  
  return (0);
}


void *EcaR75p(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_r75p", 0, 31, NULL);

  request.var1.name     = R75P_NAME;
  request.var1.longname = R75P_LONGNAME;
  request.var1.units    = R75P_UNITS;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = 1.0;
  request.var1.f2       = NULL;
  request.var1.f3       = farselgt;
  request.var1.f4       = farnum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TIME;
  request.var2.h2       = NULL;
    
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaR75ptot(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_r75ptot", 0, 31, NULL);

  request.var1.name     = R75PTOT_NAME;
  request.var1.longname = R75PTOT_LONGNAME;
  request.var1.units    = R75PTOT_UNITS;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = 1.0;
  request.var1.f2       = NULL;
  request.var1.f3       = farselgt;
  request.var1.f4       = farsum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TOTAL_AMOUNT;
  request.var2.h2       = NULL;
    
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaR90p(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_r90p", 0, 31, NULL);

  request.var1.name     = R90P_NAME;
  request.var1.longname = R90P_LONGNAME;
  request.var1.units    = R90P_UNITS;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = 1.0;
  request.var1.f2       = NULL;
  request.var1.f3       = farselgt;
  request.var1.f4       = farnum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TIME;
  request.var2.h2       = NULL;
    
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaR90ptot(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_r90ptot", 0, 31, NULL);

  request.var1.name     = R90PTOT_NAME;
  request.var1.longname = R90PTOT_LONGNAME;
  request.var1.units    = R90PTOT_UNITS;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = 1.0;
  request.var1.f2       = NULL;
  request.var1.f3       = farselgt;
  request.var1.f4       = farsum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TOTAL_AMOUNT;
  request.var2.h2       = NULL;
    
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaR95p(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_r95p", 0, 31, NULL);

  request.var1.name     = R95P_NAME;
  request.var1.longname = R95P_LONGNAME;
  request.var1.units    = R95P_UNITS;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = 1.0;
  request.var1.f2       = NULL;
  request.var1.f3       = farselgt;
  request.var1.f4       = farnum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TIME;
  request.var2.h2       = NULL;
    
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaR95ptot(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_r95ptot", 0, 31, NULL);

  request.var1.name     = R95PTOT_NAME;
  request.var1.longname = R95PTOT_LONGNAME;
  request.var1.units    = R95PTOT_UNITS;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = 1.0;
  request.var1.f2       = NULL;
  request.var1.f3       = farselgt;
  request.var1.f4       = farsum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TOTAL_AMOUNT;
  request.var2.h2       = NULL;
    
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaR99p(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_r99p", 0, 31, NULL);

  request.var1.name     = R99P_NAME;
  request.var1.longname = R99P_LONGNAME;
  request.var1.units    = R99P_UNITS;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = 1.0;
  request.var1.f2       = NULL;
  request.var1.f3       = farselgt;
  request.var1.f4       = farnum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TIME;
  request.var2.h2       = NULL;
    
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaR99ptot(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_r99ptot", 0, 31, NULL);

  request.var1.name     = R99PTOT_NAME;
  request.var1.longname = R99PTOT_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = 1.0;
  request.var1.f2       = NULL;
  request.var1.f3       = farselgt;
  request.var1.f4       = farsum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TOTAL_AMOUNT;
  request.var2.h2       = NULL;
    
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaRr1(void *argument)
{
  ECA_REQUEST_1 request;
  char lnamebuffer[1024];
  double threshold = 1;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_rr1", 0, 31, NULL);

  if ( operatorArgc() == 1 ) threshold = parameter2double(operatorArgv()[0]);
  if ( operatorArgc() > 1 ) cdoAbort("Too many arguments!");

  sprintf(lnamebuffer, RR1_LONGNAME, threshold);

  request.var1.name     = RR1_NAME;
  request.var1.longname = lnamebuffer;
  request.var1.units    = RR1_UNITS;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = threshold;
  request.var1.f2       = farnum;
  request.var1.f3       = NULL;
  request.var1.mulc     = 0.0;    
  request.var1.addc     = 0.0; 
  request.var1.epilog   = NONE;   
  request.var2.h2       = NULL;    
  request.var2.h3       = NULL;    
   
  eca1(&request);
  cdoFinish();
  
  return (0);
}


void *EcaRx1day(void *argument)
{
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  if ( operatorArgc() > 0 )
    {
      if ( 'm' == operatorArgv()[0][0] )
        cdoOperatorAdd("eca_rx1day", 0, 8,  NULL); /* monthly mode */
      else
        cdoWarning("Parameter value '%s' is invalid. The only valid value is 'm' indicating monthly mode. Operating in yearly mode now.");
    }
  else 
    cdoOperatorAdd("eca_rx1day", 0, 31, NULL);

  request.var1.name     = RX1DAY_NAME;
  request.var1.longname = RX1DAY_LONGNAME;
  request.var1.units    = RX1DAY_UNITS;
  request.var1.f1       = NULL;
  request.var1.f2       = farmax;
  request.var1.f3       = NULL;
  request.var1.mulc     = 0.0;    
  request.var1.addc     = 0.0;    
  request.var1.epilog   = NONE;
  request.var2.h2       = NULL;
  request.var2.h3       = NULL;
   
  eca1(&request);
  cdoFinish();
  
  return (0);
}


void *EcaRx5day(void *argument)
{
  char *longname;
  double argX = 50.0;
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  if ( operatorArgc() > 0 ) argX = parameter2double(operatorArgv()[0]);
  
  longname = (char*) malloc(strlen(RX5DAY_LONGNAME2) + 40);
  sprintf(longname, RX5DAY_LONGNAME2, argX);
  
  cdoOperatorAdd("eca_rx5day", 0, 31, NULL);

  request.var1.name     = RX5DAY_NAME;
  request.var1.longname = RX5DAY_LONGNAME;
  request.var1.units    = RX5DAY_UNITS;
  request.var1.f1       = NULL;
  request.var1.f2       = farmax;
  request.var1.f3       = NULL;
  request.var1.mulc     = 0.0; 
  request.var1.addc     = 0.0; 
  request.var1.epilog   = NONE;
  request.var2.name     = RX5DAY_NAME2;
  request.var2.longname = longname;
  request.var2.units    = RX5DAY_UNITS2;
  request.var2.h1       = farselgec;
  request.var2.h1arg    = argX;
  request.var2.h2       = farnum;
  request.var2.h3       = NULL;
   
  eca1(&request);
  
  free(longname);
  cdoFinish();
  
  return (0);
}


void *EcaSdii(void *argument)
{
  ECA_REQUEST_1 request;
  char lnamebuffer[1024];
  double threshold = 1;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_sdii", 0, 31, NULL);

  if ( operatorArgc() == 1 ) threshold = parameter2double(operatorArgv()[0]);
  if ( operatorArgc() > 1 ) cdoAbort("Too many arguments!");

  sprintf(lnamebuffer, SDII_LONGNAME, threshold);

  request.var1.name     = SDII_NAME;
  request.var1.longname = lnamebuffer;
  request.var1.units    = SDII_UNITS;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = threshold;
  request.var1.f2       = farsum;
  request.var1.f3       = NULL;
  request.var1.mulc     = 0.0;    
  request.var1.addc     = 0.0;    
  request.var1.epilog   = MEAN;
  request.var2.h2       = NULL;
  request.var2.h3       = NULL;
   
  eca1(&request);
  cdoFinish();
  
  return (0);
}


void *Fdns(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("fdns", 0, 31, NULL);

  request.var1.name     = FDNS_NAME;
  request.var1.longname = FDNS_LONGNAME;
  request.var1.units    = FDNS_UNITS;
  request.var1.f1       = farsellec;
  request.var1.f1arg    = TO_KELVIN(0.0);
  request.var1.f2       = farsellec;
  request.var1.f2arg    = 0.01;
  request.var1.f3       = faradd; /* any f with f(a, b) = miss, if a = miss or b = miss will do here */
  request.var1.f4       = farnum;
  request.var1.f5       = NULL;
  request.var1.epilog   = NONE;
  request.var2.h2       = NULL;
    
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *Strwin(void *argument)
{
  char *longname;
  double maxWind = 10.5;
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("strwin", 0, 31, NULL);

  if ( operatorArgc() > 0 )
    maxWind = parameter2double(operatorArgv()[0]);

  longname = (char*) malloc(strlen(STRWIN_LONGNAME) + 40);
  sprintf(longname, STRWIN_LONGNAME, maxWind);
         
  request.var1.name     = STRWIN_NAME;
  request.var1.longname = longname;
  request.var1.units    = STRWIN_UNITS;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = maxWind;
  request.var1.f2       = farnum;
  request.var1.f3       = NULL;
  request.var1.mulc     = 0.0;    
  request.var1.addc     = 0.0;    
  request.var1.epilog   = NONE;
  request.var2.name     = STRWIN_NAME2;
  request.var2.longname = STRWIN_LONGNAME2;
  request.var2.units    = STRWIN_UNITS2;
  request.var2.h1       = farselgec;
  request.var2.h1arg    = maxWind;
  request.var2.h2       = farnum2;
  request.var2.h3       = farmax;
   
  eca1(&request);
  
  free(longname);
  cdoFinish();
  
  return (0);
}


void *Strbre(void *argument)
{
  static const double maxWind = 10.5;
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("strbre", 0, 31, NULL);
         
  request.var1.name     = STRBRE_NAME;
  request.var1.longname = STRBRE_LONGNAME;
  request.var1.units    = STRWIN_UNITS;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = maxWind;
  request.var1.f2       = farnum;
  request.var1.f3       = NULL;
  request.var1.mulc     = 0.0;    
  request.var1.addc     = 0.0;    
  request.var1.epilog   = NONE;
  request.var2.name     = STRBRE_NAME2;
  request.var2.longname = STRBRE_LONGNAME2;
  request.var2.units    = STRWIN_UNITS2;
  request.var2.h1       = farselgec;
  request.var2.h1arg    = maxWind;
  request.var2.h2       = farnum2;
  request.var2.h3       = farmax;
   
  eca1(&request);
  cdoFinish();
  
  return (0);
}


void *Strgal(void *argument)
{
  static const double maxWind = 20.5;
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("strgal", 0, 31, NULL);
         
  request.var1.name     = STRBRE_NAME;
  request.var1.longname = STRBRE_LONGNAME;
  request.var1.units    = STRWIN_UNITS;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = maxWind;
  request.var1.f2       = farnum;
  request.var1.f3       = NULL;
  request.var1.mulc     = 0.0;    
  request.var1.addc     = 0.0;    
  request.var1.epilog   = NONE;
  request.var2.name     = STRBRE_NAME2;
  request.var2.longname = STRBRE_LONGNAME2;
  request.var2.units    = STRWIN_UNITS2;
  request.var2.h1       = farselgec;
  request.var2.h1arg    = maxWind;
  request.var2.h2       = farnum2;
  request.var2.h3       = farmax;
   
  eca1(&request);
  cdoFinish();
  
  return (0);
}


void *Hurr(void *argument)
{
  static const double maxWind = 32.5;
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("hurr", 0, 31, NULL);
         
  request.var1.name     = HURR_NAME;
  request.var1.longname = HURR_LONGNAME;
  request.var1.units    = STRWIN_UNITS;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = maxWind;
  request.var1.f2       = farnum;
  request.var1.f3       = NULL;
  request.var1.mulc     = 0.0;    
  request.var1.addc     = 0.0;    
  request.var1.epilog   = NONE;
  request.var2.name     = HURR_NAME2;
  request.var2.longname = HURR_LONGNAME2;
  request.var2.units    = STRWIN_UNITS2;
  request.var2.h1       = farselgec;
  request.var2.h1arg    = maxWind;
  request.var2.h2       = farnum2;
  request.var2.h3       = farmax;
   
  eca1(&request);
  cdoFinish();
  
  return (0);
}
