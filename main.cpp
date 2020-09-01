#include "Membrane_Helfrich.hpp"
double L=1000,radius;
double MOVER=1.0;
double K=140.0,H_0=0.02;
double K_G=0.0;
double Ka=4000000,Surface_A0=350000;
double Kv=0.0,eta=1.0;

//(L,radius,MOVER,K,H_0,K_G,Ka,dP);

int main() {
    //int WS=1000000;
    //int MS=2000000;
    int WS=5000000,MS=1000000;//For debug
    //double Factor=1.2;
    //double Factor0=10;

    //int N_Annual=50;
    //int Cycles_Per_Sample= 100000;
    int Cycles_Per_Sample= 10000;
    int SampleT=0;

	ifstream TestIn("Tri.txt");
	int2 temp;
	int3 Ltemp;
	temp.y=-1;
	int Ftest=0;
	while(!TestIn.eof()) {
		TestIn>>temp.x>>Ltemp.x>>Ltemp.y>>Ltemp.z;
		if(temp.y==temp.x)break;
		temp.y=temp.x;
		Ftest++;
	}
	TestIn.close();
	cout<<"F is"<<Ftest<<endl;
	radius=sqrt(Surface_A0/(1.7*sqrt(3)*Ftest));
	//K*=3;
	//Ka*=3;
	Membrane_Helfrich Membrane_Simu(L, radius, MOVER, K, H_0, K_G, Ka, Kv, Surface_A0, eta);

    int F=Membrane_Simu.getF();
    int V=Membrane_Simu.getV();
    int E=Membrane_Simu.getE();


    for(int ms=0;ms<(WS);ms++) {
        double acc;
        bool flag=false;

        acc=0;
        for(int l=0;l<V;l++)acc+=Membrane_Simu.move_node(l);
        acc/=V;
        if((acc<0.4)||(acc>0.6))cout<<"acc(move node)"<<acc<<endl;
        if(acc<0.4){MOVER*=(acc/0.5);flag=true;}
        if(acc>0.6){MOVER*=(acc/0.5);flag=true;}
        
        if(flag)Membrane_Simu.SetTrialPara(MOVER);


        if(ms%200==0){
                acc=0;
                for(int l=0;l<E;l++)acc+=Membrane_Simu.flip_edge(l);
                acc/=E;
                cout<<"acc(flip)"<<acc<<endl;
                if(ms%10000==0)Membrane_Simu.Check();//FOR DEBUG
                Membrane_Simu.Refresh_Sys_Var();
                if(ms%10000==0)cout<<"Warming "<<ms*100.0/(WS)<<"%"<<endl;
                if(ms%10000==0)Membrane_Simu.Print(0);//For test.
        }
    }

//    K*=3;
//    Ka*=3;
    Membrane_Simu.SetMembranePara(K,K_G,Ka,Kv);

    for(int ms=0;ms<(WS);ms++) {
        double acc;
        bool flag=false;

        acc=0;
        for(int l=0;l<V;l++)acc+=Membrane_Simu.move_node(l);
        acc/=V;
        if((acc<0.4)||(acc>0.6))cout<<"acc(move node)"<<acc<<endl;
        if(acc<0.4){MOVER*=(acc/0.5);flag=true;}
        if(acc>0.6){MOVER*=(acc/0.5);flag=true;}
        
        if(flag)Membrane_Simu.SetTrialPara(MOVER);


        if(ms%200==0){
                acc=0;
                for(int l=0;l<E;l++)acc+=Membrane_Simu.flip_edge(l);
                acc/=E;
                cout<<"acc(flip)"<<acc<<endl;
                if(ms%10000==0)Membrane_Simu.Check();//FOR DEBUG
                Membrane_Simu.Refresh_Sys_Var();
                if(ms%10000==0)cout<<"Warming "<<ms*100.0/(WS)<<"%"<<endl;
                if(ms%10000==0)Membrane_Simu.Print(1);//For test.
        }
    }


//    K*=Factor0;
//    K_G*=Factor0;
//    Ka*=Factor0;
//    Kv*=Factor0;
//    Membrane_Simu.SetMembranePara(K,K_G,Ka,Kv);

/*
    for(int I_An=0;I_An<N_Annual;I_An++) {
        for(int ms=0;ms<WS;ms++) {
            double acc;
            bool flag=false;

            acc=0;
            for(int l=0;l<V;l++)acc+=Membrane_Simu.move_node(l);
            acc/=V;
            if((acc<0.3)||(acc>0.6))cout<<"acc(move node)"<<acc<<endl;
            if(acc<0.3){MOVER*=(acc/0.5);flag=true;}
            if(acc>0.6){MOVER*=(acc/0.5);flag=true;}
        
            if(flag)Membrane_Simu.SetTrialPara(MOVER);


            if(ms%200==0){
                acc=0;
                for(int l=0;l<E;l++)acc+=Membrane_Simu.flip_edge(l);
                acc/=E;
                cout<<"acc(flip)"<<acc<<endl;
            	if(ms%10000==0)Membrane_Simu.Check();//FOR DEBUG
            	Membrane_Simu.Refresh_Sys_Var();
                if(ms%10000==0)cout<<"Warming "<<ms*100.0/WS<<"%"<<endl;
            }
        }
    

        for(int ms=0;ms<MS;ms++) {
            for(int l=0;l<V;l++)Membrane_Simu.move_node(l);

            if(ms%200==0){
                for(int l=0;l<E;l++)Membrane_Simu.flip_edge(l);
                //if(ms%10000==0)Membrane_Simu.Check();//FOR DEBUG
            	Membrane_Simu.Refresh_Sys_Var();
            }
            if(ms%Cycles_Per_Sample==0){
                SampleT++;
                Membrane_Simu.Print(SampleT);cout<<"Round "<<I_An<<" "<<ms<<"cycles"<<endl;
            }
        }

        K*=Factor;
        K_G*=Factor;
        Ka*=Factor;
        Kv*=Factor;
        //Membrane_Simu.Check();
        Membrane_Simu.SetMembranePara(K,K_G,Ka,Kv);

    }*/
}
