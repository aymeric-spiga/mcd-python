function DefaultDateValues(){

var rightnow= new Date();

//document.calendar.year.value=2000;
document.calendar.year.value=rightnow.getFullYear();
//document.calendar.month.value=1;
document.calendar.month.value=rightnow.getMonth()+1;
//document.calendar.day.value=1;
document.calendar.day.value=rightnow.getDate();
}

function DefaultTimeValues(){

var rightnow= new Date();

document.calendar.hours.value=rightnow.getUTCHours();
//document.calendar.hours.value=0;
document.calendar.minutes.value=rightnow.getUTCMinutes();
//document.calendar.minutes.value=0;
document.calendar.seconds.value=rightnow.getUTCSeconds();
//document.calendar.seconds.value=0;
}

function DateValues(y,m,d){
var y;
var m;
var d;
document.calendar.year.value=y;
document.calendar.month.value=m;
document.calendar.day.value=d;
}

function TimeValues(h,m,s){
var h;
var m;
var s;
document.calendar.hours.value=h;
document.calendar.minutes.value=m;
document.calendar.seconds.value=s;
}

function DateAndTimeValues(year,month,day,hours,minutes,seconds){
var year,month,day,hours,minutes,seconds;
DateValues(year,month,day);
TimeValues(hours,minutes,seconds);
}

function CheckGivenYear(){
var bissextil; // bissextil year ? (0==no, 1==yes) (returned value)
var val=document.calendar.year.value;

while ((val!=Math.round(val))||(val<1583)) {
  if (val!=Math.round(val)) {
    val=prompt("Year must be an integer! e.g.:",Math.round(val));
  }
  if (val<1583) {
    val=prompt("Sorry, Year must not be less than 1583 \n (in order to stick to the Gregorian calendar)",1583);
  }
}

document.calendar.year.value=val;

// check if it is a bissextil year
/* a year is bissextil if it is a multiple of 4 but not of 100,
or if it is a multiple of 400 */
if ((((val%4)==0)&&((val%100)!=0))||((val%400)==0)) {
  bissextil=1;
}
else {
  bissextil=0; // not a bissextil year
}

return bissextil;
}

function CheckGivenMonth(){
var val=document.calendar.month.value;

while((val<=0)||(val>12)||(val!=Math.round(val))) {
 val=prompt("Month must be an integer between 1 and 12!",Math.round(val));
}

document.calendar.month.value=val;
}

function CheckGivenDay(year,bissextil,month){
var year;
var bissextil; // bissextil year ? (0==no, 1==yes)
var month;
var dayval;

dayval=document.calendar.day.value;

if (month==1) {
  while (((dayval<=0)||(dayval>31))||(dayval!=Math.round(dayval))){
    dayval=prompt("Invalid day! (must be an integer ranging from 1 to 31)",dayval);
  }
}
if (month==2) {
  if (bissextil==0){
    while (((dayval<=0)||(dayval>28))||(dayval!=Math.round(dayval))){
      dayval=prompt("Invalid day! (must be an integer ranging from 1 to 28)",dayval);
    }
  }
  else {
    while (((dayval<=0)||(dayval>29))||(dayval!=Math.round(dayval))){
      dayval=prompt("Invalid day! (must be an integer ranging from 1 to 29)",dayval);
    }
  }
}
if (month==3) {
  while (((dayval<=0)||(dayval>31))||(dayval!=Math.round(dayval))){
    dayval=prompt("Invalid day! (must be an integer ranging from 1 to 31)",dayval);
  }
}
if (month==4) {
  while (((dayval<=0)||(dayval>30))||(dayval!=Math.round(dayval))){
    dayval=prompt("Invalid day! (must be an integer ranging from 1 to 30)",dayval);
  }
}
if (month==5) {
  while (((dayval<=0)||(dayval>31))||(dayval!=Math.round(dayval))){
    dayval=prompt("Invalid day! (must be an integer ranging from 1 to 31)",dayval);
  }
}
if (month==6) {
  while (((dayval<=0)||(dayval>30))||(dayval!=Math.round(dayval))){
    dayval=prompt("Invalid day! (must be an integer ranging from 1 to 30)",dayval);
  }
}
if (month==7) {
  while (((dayval<=0)||(dayval>31))||(dayval!=Math.round(dayval))){
    dayval=prompt("Invalid day! (must be an integer ranging from 1 to 31)",dayval);
  }
}
if (month==8) {
  while (((dayval<=0)||(dayval>31))||(dayval!=Math.round(dayval))){
    dayval=prompt("Invalid day! (must be an integer ranging from 1 to 31)",dayval);
  }
}
if (month==9) {
  while (((dayval<=0)||(dayval>30))||(dayval!=Math.round(dayval))){
    dayval=prompt("Invalid day! (must be an integer ranging from 1 to 30)",dayval);
  }
}
if (month==10) {
  while (((dayval<=0)||(dayval>31))||(dayval!=Math.round(dayval))){
    dayval=prompt("Invalid day! (must be an integer ranging from 1 to 31)",dayval);
  }
}
if (month==11) {
  while (((dayval<=0)||(dayval>30))||(dayval!=Math.round(dayval))){
    dayval=prompt("Invalid day! (must be an integer ranging from 1 to 30)",dayval);
  }
}
if (month==12) {
  while (((dayval<=0)||(dayval>31))||(dayval!=Math.round(dayval))){
    dayval=prompt("Invalid day! (must be an integer ranging from 1 to 31)",dayval);
  }
}


document.calendar.day.value=dayval;
}

function CheckGivenTime(){
var hours;
var minutes;
var seconds;

// Check value of hours
hours=document.calendar.hours.value;
while (((hours<0)||(hours>=24))||(hours!=Math.round(hours))) {
  hours=prompt("Invalid Time (hh) value! (must be an integer ranging from 0 to 23)",hours);
  document.calendar.hours.value=hours;
}

// Check value of minutes
minutes=document.calendar.minutes.value;
while (((minutes<0)||(hours>=60))||(minutes!=Math.round(minutes))) {
  minutes=prompt("Invalid Time (mm) value! (must be an integer ranging from 0 to 59)",minutes);
  document.calendar.minutes.value=minutes;
}

//Check value of seconds
seconds=document.calendar.seconds.value;
while((seconds<0)||(seconds>=60)) {
  seconds=prompt("Invalid Time (ss) value! (must range from 0 to 60)",seconds);
  document.calendar.seconds.value=seconds;
}

}

/*--------------------------------------------------------------------------*/

function CheckGivenDate(){
var bissextil; // bissextil year ? (0==no, 1==yes)

bissextil=CheckGivenYear();
CheckGivenMonth();
CheckGivenDay(document.calendar.year.value,bissextil,document.calendar.month.value);
CheckGivenTime();
//alert("OK");
return bissextil;
}

/*--------------------------------------------------------------------------*/

function Convert2Julian(){
var bissextil; // bissextil year ? (0==no, 1==yes)
var year;
var month;
var day;
var i;
var hours,minutes,seconds;
var ref_year=1968;
var ref_jdate=2.4398565e6; // Julian date for 01/01/1968 00:00:00
var edays = new Array(0,31,59,90,120,151,181,212,243,273,304,334);
// edays = number of elapsed days during previous monthes of same year
var nday=0.0; // number of days

// start by checking validity of given date
bissextil=CheckGivenDate();

year=document.calendar.year.value;
month=document.calendar.month.value;
day=document.calendar.day.value;

// compute number of days due to years 
if(year>ref_year) {
  for(i=ref_year;i<year;i++){
    nday=nday+365.0;
    if ((((i%4)==0)&&((i%100)!=0))||((i%400)==0)) { // bissextil year
      nday++;
    }
  }
}
else {
  for(i=year;i<ref_year;i++){
    nday=nday-365.0;
    if ((((i%4)==0)&&((i%100)!=0))||((i%400)==0)) { // bissextil year
      nday--;
    }
  }
}

// add number of days due to elapsed monthes
nday=nday+edays[month-1];
//alert(nday)

//add 1 if year is bissextil and month >=3
if((bissextil==1)&&(month>=3)){
  nday=nday+1;
}
// add reference year offset and day
//jdate=ref_jdate+nday+day;
jdate=nday*1.0+day*1.0+ref_jdate*1.0-1.0;

// add time (hours+minutes+seconds)
hours=document.calendar.hours.value;
minutes=document.calendar.minutes.value;
seconds=document.calendar.seconds.value;
jdate=jdate+hours/24.0+minutes/1440.0+seconds/86400.0;

document.calendar.julian.value=jdate;
}

/*--------------------------------------------------------------------------*/

function Convert2Ls(){
// Convert a Julian date to corresponding "sol" and "Ls"
var jdate;
var sol;
var ls;
var martianyear;
var martianmonth;

var jdate_ref=2.442765667e6; // 19/12/1975 4:00:00, such that Ls=0
// jdate_ref is also the begining of Martian Year "12"
var martianyear_ref=12;
var earthday=86400.0;
var marsday=88775.245;
var marsyear=668.60; // number of sols in a martian year 

// Start by converting given date to Julian date
Convert2Julian();

// Convert julian days to sol date
jdate=document.calendar.julian.value;

sol=(jdate-jdate_ref)*earthday/marsday;

martianyear=martianyear_ref;
// Compute Martian Year #, along with sol value
// sol being computed modulo the number of sols in a martian year
while (sol>=marsyear){
  sol=sol-marsyear;
  martianyear=martianyear+1;
}
while (sol<0.0){
  sol=sol+marsyear;
  martianyear=martianyear-1;
}

//document.dummy.dummy1.value=sol;

// convert sol number to Ls
ls=Sol2Ls(sol);

// Knowing Ls compute martian month
martianmonth=1+Math.floor(ls/30.);

//Display value with a maximum of 2 decimal digits
document.calendar.martianyear.value=martianyear;
document.calendar.martianmonth.value=martianmonth;
document.calendar.ls.value=Math.round(ls*10)/10;
//document.calendar.sol.value=Math.round(sol*10)/10;
document.calendar.sol.value=1+Math.floor(sol);
}

/*--------------------------------------------------------------------------*/

function Sol2Ls(sol) {
var sol;
var ls;

var year_day=668.6; // number of sols in a martian year
var peri_day=485.35; // perihelion date
var e_ellip=0.09340; // orbital ecentricity
var timeperi=1.90258341759902 // 2*Pi*(1-Ls(perihelion)/360); Ls(perihelion)=250.99
var rad2deg=180./Math.PI;

var i;
var zz,zanom,zdx=10;
var xref,zx0,zteta;
// xref: mean anomaly, zx0: eccentric anomaly, zteta: true anomaly

zz=(sol-peri_day)/year_day;
zanom=2.*Math.PI*(zz-Math.round(zz));
xref=Math.abs(zanom);

// Solve Kepler equation zx0 - e *sin(zx0) = xref
// Using Newton iterations
zx0=xref+e_ellip*Math.sin(xref);
do {
  zdx=-(zx0-e_ellip*Math.sin(zx0)-xref)/(1.-e_ellip*Math.cos(zx0));
  zx0=zx0+zdx;
}while (zdx>1.e-7);
if (zanom<0) zx0=-zx0;

// Compute true anomaly zteta, now that eccentric anomaly zx0 is known
zteta=2.*Math.atan(Math.sqrt((1.+e_ellip)/(1.-e_ellip))*Math.tan(zx0/2.));

// compute Ls
ls=zteta-timeperi;
if(ls<0) ls=ls+2.*Math.PI;
if(ls>2.*Math.PI) ls=ls-2.*Math.PI;
// convert Ls into degrees
ls=rad2deg*ls;

return ls;
}


function CorrectValue(){
if ( document.calendar.longitude.value == "all") {
  document.calendar.longitude.value = 0.
}
if ( document.calendar.latitude.value == "all") {
  document.calendar.latitude.value = 0.
}
}

function PlaceValues(lon,lat){
var lon;
var lat;
var savlon;
var savlat;
var yeah;
savlon=document.calendar.longitude.value;
savlat=document.calendar.latitude.value;
document.calendar.longitude.value=lon;
document.calendar.latitude.value=lat;
yeah=document.calendar.toto.value;
// to be improved. lon and lat should be separated. but for the moment
// only mapping capabilities exist as 2D plots which means if lon is all
// lat is all and vice versa.
if ( document.calendar.longitude.value == "all") {
 if ( document.calendar.altitude.value == "all") {
    alert("Not allowed! Altitude and time must be fixed in mapping mode.");
    document.calendar.longitude.value=savlon;
    document.calendar.latitude.value=savlat;
    document.calendar.toto.value=yeah;
 }
 if ( document.calendar.localtime.value == "all") {
    alert("Not allowed! Altitude and time must be fixed in mapping mode.");
    document.calendar.longitude.value=savlon;
    document.calendar.latitude.value=savlat;
    document.calendar.toto.value=yeah;
 }
}
}

function PlaceVar(var1,var2,var3,var4){
var var1;
var var2;
var var3;
var var4;
document.calendar.var1.value=var1;
document.calendar.var2.value=var2;
document.calendar.var3.value=var3;
document.calendar.var4.value=var4;
}

function DefaultSpaceTime(){
document.calendar.localtime.value = 9
document.calendar.altitude.value = 2
}

// this is separated from PlaceValues because
// allows easier determination through earth time
// if needed (e.g. example case with rovers)
function PlaceValues2(localtime){
var localtime;
var savlocaltime;
savlocaltime=document.calendar.localtime.value
document.calendar.localtime.value=localtime;
//if ( document.calendar.localtime.value == "all" ) {
// if ( document.calendar.altitude.value == "all" ) {
//    alert("Not allowed! Now choose either a fixed altitude or a fixed time.");
//    //document.calendar.localtime.value=savlocaltime;
// }
//}
//if ( document.calendar.localtime.value == "all" ) {
// if ( document.calendar.longitude.value == "all" ) {
//    alert("Not allowed! Now choose either a lander or a fixed time.");
//    //document.calendar.localtime.value=savlocaltime;
// }
//}
}
function PlaceValues3(altitude){
var altitude;
var savaltitude;
savaltitude=document.calendar.altitude.value
document.calendar.altitude.value=altitude;
//if ( document.calendar.altitude.value == "all" ) {
// if ( document.calendar.localtime.value == "all") {
//    alert("Not allowed! Now choose either a fixed time or a fixed altitude.");
//    //document.calendar.altitude.value=savaltitude;
// }
//}
//if ( document.calendar.altitude.value == "all" ) {
// if ( document.calendar.longitude.value == "all") {
//    alert("Not allowed! Now choose either a lander or a fixed altitude.");
//    //document.calendar.altitude.value=savaltitude;
// }
//}
}

function Checkaltitude(){
altitude=document.calendar.altitude.value
num=+altitude // returns numeric value, or NaN if string not purely numeric characters
while ( isNaN(num) ) {
 num=prompt("\nAltitude is set as 'all' or a range, which makes it a free dimension.\nPlease indicate a numerical value to make the requested mode work:\n",1000)
}
document.calendar.altitude.value=num
}

function Checklocaltime(){
localtime=document.calendar.localtime.value
num=+localtime // returns numeric value, or NaN if string not purely numeric characters
while ( isNaN(num) ) {
 num=prompt("\nLocal time is set as 'all' or a range, which makes it a free dimension.\nPlease indicate a numerical value to make the requested mode work:\n",12)
}
document.calendar.localtime.value=num
}

function Checklongitude(){
longitude=document.calendar.longitude.value
num=+longitude // returns numeric value, or NaN if string not purely numeric characters
while ( isNaN(num) ) {
 num=prompt("\nLongitude is set as 'all' or a range, which makes it a free dimension.\nPlease indicate a numerical value to make the requested mode work:\n",0)
}
document.calendar.longitude.value=num
}

function Checklatitude(){
latitude=document.calendar.latitude.value
num=+latitude // returns numeric value, or NaN if string not purely numeric characters
while ( isNaN(num) ) {
 num=prompt("\nLatitude is set as 'all' or a range, which makes it a free dimension.\nPlease indicate a numerical value to make the requested mode work:\n",0)
}
document.calendar.latitude.value=num
}


function ChooseColor(colorm){
var colorm;
document.calendar.colorm.value=colorm;
}

function submit_form_beginner() {
submit_form()
document.calendar.reset();
DefaultDateValues();
DefaultTimeValues();
Convert2Ls();
PlaceValues(0.,0.);
DefaultSpaceTime();
}

function submit_form() {
Convert2Julian();
document.calendar.submit();
}





