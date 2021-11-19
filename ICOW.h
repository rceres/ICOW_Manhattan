

void iCOWmodule(double* levers, double* objs, double* info, double* consts);
void CharacterizeCity (double W,double B,double R,double P, double D, double * cityChar);


const double costUnit=1000000000;
const double CEC=17;              //m City Elevation Change, Bennet Park in the Washington Heights area of Manhattan is 
const double CityWidth=43000.0;                      //m
const double CityLength=2000.0;                     //m
const double TotalCityValueInitial = 1500000000000; // 1,500,000,000,000 1,000,000,000,000  1,00,000,000,000;  50,000,000,000,000
const double BH = 24;//30;//20; //m

const double CitySlope=CityLength/CityWidth;
const double WithdrawelPercentLost = 0.01;
const double ProtectedValueRatio = 1.1;
const double SlopeDike = .5;
const double DikeUnprotectedValuationRatio = 1.0;
const double WidthDikeTop = 3; //m
const double DikeStartingCostPoint = 2;
const double UnitCostPerVolumeDike = 10; //$ dollars per m^3

const double WithdrawelCostFactor = 1.0;
const double resistanceExponentialFactor = 0.115;
const double resistanceLinearFactor=0.35;
const double resistanceExponentialThreshold = .4;
const double damageFactor = 0.39;
// i.e. damage is worse when the dike fails
const double FailedDikeDamageFactor = 1.5; // considers additional damage that resutls because of dike failure
const double intactDikeDamageFactor = 0.03;
const double pfThreshold=0.95;
const double pfBase=.05;
const double minHeight=0.1;
const double WminHeight=0.01;
const double Basement=3.0;
const double threshold = TotalCityValueInitial/375;
const double thresholdDamageFraction = 0.5;
// threhold is a demarcation of damage that is considered unacceptable
// thresholdDamageFraction = 0 causes damage to accumulate at the normal (below threshold) rate
// threshioldDamageMultiple = 1 causes damage to accululate at normal + normal (2x) below threshold rate
const double thresholdDamageExponent = 1.01;

const int lengthSurgeSequences=200;
   
const double baseValue=100;
const double PBase=0.5;
const double Seawall=2;  // from Talke
const double runUpWave=1.1; // to account for wave/runup. 1.0 results in no increase
const int maxSurgeBlock=5000;


// index lables for the city
const int caseNum=0;
const int wh=1;
const int rh=2;
const int rp=3;
const int dbh=4;
const int dh=5;
const int vz1=6;
const int vz2=7;
const int vz3=8;
const int vz4=9;
const int tz1=10;
const int tz2=11;
const int tz3=12;
const int tz4=13;
const int fw=14;
const int tcvi=15; // total city value initial
const int ilfw=16;
const int tcvaw=17;
const int vifod=18;
const int vbd=19;
const int fcv=20;
const int dc=21;
const int wc=22;
const int rc=23;
const int tic=24; // total investment cost
const int tc=25; // total net cost includes TIC plus loss of city value
const int dtr=26;
const int numCityChar=27;

// index values for the damageVector
const int dvt=0;  // total damage cost
const int dvz1=1;  // damage zone 1
const int dvz2=2;  // damage zone 2
const int dvz3=3;  // damage zone 3
const int dvz4=4;  // damage zone 4
const int dvFE=5;  // Flood Event, some damage occurs
const int dvBE=6;  // Breech Event
const int dvTE=7;  // Threshold Event
const int dvLength=8;

const double PL = 0.05;
const double PH = 0.995;
const double DHL=0;
const double DHH=17;
const double RHL=0;
const double RHH=14;
const double WHL=0;
const double WHH=6;
const double DBL=0;
const double DBH=6;
const double WHD=0;

const int years=50;
int statesToEvaluate=5000;
const int maxStatesToEvaluate=20000;
double surges [lengthSurgeSequences*maxStatesToEvaluate]; // years*stateToEvaluate

unsigned m_w=521288629;
unsigned m_z=362436069;
double maxBaseSurge=12;
double surgeMaxRateIncrease=0;
