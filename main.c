#include "ellipsoid_fitting.h"
#include <stdio.h>
double data1[]={-147,-96,179,-157,-89,198,-153,-98,194,-168,-85,163,-165,-80,165,-168,-80,170,-169,-71,179,-181,-58,179,-180,-64,165,-180,-64,105,-136,-75,25,-72,-72,-32,17,-50,-43,144,-28,24,176,58,19,155,108,-18,47,131,-76,-127,157,-48,-201,180,67,-166,249,125,-85,295,174,-51,285,239,57,220,306,191,117,251,228,63,116,229,44,119,118,-14,337,-4,-53,351,-68,-66,341,-202,-52,219,-180,-62,20,-67,-95,-39,89,-70,-27,186,-58,92,137,-73,242,91,-80,299,15,-67,337,-89,-54,328,-126,-70,308,6,-54,347,121,-44,284,171,-73,72,111,-92,1,11,-98,-20,-83,-119,-4,-196,-84,167,-191,-66,239,-203,5,251,-245,109,179,-180,189,18,-20,309,49,-13,319,190,-126,212,299,-133,59,349,-76,-7,360,137,-66,280,174,-38,55,91,-72,-30,195,-62,134,116,-67,308,-125,-35,329,-191,-43,249,-152,-99,41,-65,-95,-30,52,-90,-47,152,-96,108,170,-80,197,161,-45,263,118,8,325,105,86,346,104,118,342,89,210,325,67,238,299,51,308,196,6,322,131,7,275,2,22,241,-43,11,205,-64,22,76,-106,65,-11,-76,145,-92,32,145,-87,73,118,-82,3,143,-54,8,-41,-53,-75,-252,30,125,-233,84,227,-141,148,342,-2,152,375,183,206,210,188,163,53,159,124,-36,125,76,-67,19,38,-113,-121,8,-71,-206,-15,259,-179,-40,283,-171,-30,288,-176,-94,230,-40,-161,192,8,-164,71,89,-156,144,91,-159,168,58,-140,263,-34,-137,261,-93,-165,159,-104,-164,153,-85,-160,164,-76,-168,147,8,-166,117,73,-160,159,67,-161,105,-196,-99,121,-192,-44,260,-152,-66,278,41,-56,353,191,-48,179,133,-49,-3,44,-114,-16,96,-47,320,107,43,346,100,150,341,79,311,163,24,305,70,-17,305,50,-19,290,16,-22,101,-104,12,22,-96,-36,-71,-50,-120,-129,53,-72,-141,47,-12,-168,84,17,-157,60,-24,-89,-48,8,-44,-76,24,-28,-87,36,7,-96,25,44,-112,-20,85,-110,-50,113,-102,-49,140,-92,-34,221,-62,-22,247,-41,-26,259,-26,36,315,87,39,318,157,65,306,170,66,290,201,82,278,228};

int main(void) {
    struct Matrix data = {.entry = data1, 149, 3};
    double* dataD = malloc(sizeof(double)*200*9);
    struct Matrix D = {.entry = dataD, 149, 9};
    struct ellipsoid want;
    ellipsoid_fitting(data, D, want);
    return 0;
}