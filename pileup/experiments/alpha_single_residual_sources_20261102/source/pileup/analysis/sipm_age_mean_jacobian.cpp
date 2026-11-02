// Same six-step cross-talk closure as sipm_age_mean.cpp, with exact forward
// sensitivities of that discrete algorithm. No fast-math or physics changes.
#include <cmath>
#include <algorithm>
#include <vector>
extern "C" int gagg_age_gradient(const double* input,const double* input_gradient,
    int n,int dimensions,double cells,const double* recovery,const double* pde,
    const double* trigger,int ages,double xt,double* output,double* gradient) {
    if(n<1||dimensions<1||dimensions>8||ages<2||cells<=0)return 1;
    std::vector<double> state(ages,0.),dstate(dimensions*ages,0.),fire(ages),
        dfire(dimensions*ages),survived(ages),dsurvived(dimensions*ages),
        dcross(dimensions),dq(dimensions);
    state.back()=1.;
    for(int i=0;i<n;++i){
        if(!std::isfinite(input[i])||input[i]<0)return 2;
        double cross=0.,charge=0.;
        std::fill(dcross.begin(),dcross.end(),0.);
        for(int iteration=0;iteration<6;++iteration){
            charge=0.;std::fill(dq.begin(),dq.end(),0.);
            for(int j=0;j<ages;++j){
                double h=input[i]/cells*pde[j]+cross*trigger[j];
                fire[j]=-std::expm1(-h);
                charge+=state[j]*(fire[j]*recovery[j]);
                for(int k=0;k<dimensions;++k){
                    int a=k*ages+j;
                    dfire[a]=(1.-fire[j])*(input_gradient[k*n+i]/cells*pde[j]+dcross[k]*trigger[j]);
                    dq[k]+=(dstate[a]*fire[j]+state[j]*dfire[a])*recovery[j];
                }
            }
            cross=xt*charge;
            for(int k=0;k<dimensions;++k)dcross[k]=xt*dq[k];
        }
        output[i]=cells*charge;
        for(int k=0;k<dimensions;++k)gradient[k*n+i]=cells*dq[k];
        double fired=0.;std::fill(dq.begin(),dq.end(),0.);
        for(int j=0;j<ages;++j){
            double f=state[j]*fire[j];fired+=f;survived[j]=state[j]-f;
            for(int k=0;k<dimensions;++k){
                int a=k*ages+j;double df=dstate[a]*fire[j]+state[j]*dfire[a];
                dq[k]+=df;dsurvived[a]=dstate[a]-df;
            }
        }
        for(int j=1;j<ages;++j)state[j]=survived[j-1];
        state.back()+=survived.back();state.front()=fired;
        for(int k=0;k<dimensions;++k){
            for(int j=1;j<ages;++j)dstate[k*ages+j]=dsurvived[k*ages+j-1];
            dstate[(k+1)*ages-1]+=dsurvived[(k+1)*ages-1];dstate[k*ages]=dq[k];
        }
    }
    return 0;
}
