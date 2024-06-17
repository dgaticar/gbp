#ifndef INCLUDED_SDSL_PEMB
#define INCLUDED_SDSL_PEMB

#include <bits/stdc++.h>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/sdsl_concepts.hpp>
#include <sdsl/vectors.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v.hpp>
#include <sdsl/select_support_mcl.hpp>
#include <sdsl/wt_helper.hpp>
#include <sdsl/bp_support_sada.hpp>
#include <sdsl/inv_perm_support.hpp>
//! Namespace for the succinct data structure library.
using namespace std;
namespace sdsl
{



  template<class t_bitvector   = bit_vector,
  		class n_bitvector = bit_vector,
	   class t_succ_tree    = bp_support_sada<>,
	   class t_rank        = typename t_bitvector::rank_1_type,
	   class t_select1     = typename t_bitvector::select_1_type,
	   class t_select0     = typename t_bitvector::select_0_type,
	   class n_rank        = typename n_bitvector::rank_1_type,
	   class n_select1     = typename n_bitvector::select_1_type,
	   class n_select0     = typename n_bitvector::select_0_type>




class gbp{
    public:
        typedef int_vector<>::size_type              size_type;
        typedef int_vector<>::value_type             value_type;
        typedef random_access_const_iterator<gbp> const_iterator;
        typedef const_iterator                       iterator;
        typedef t_bitvector                          bit_vector_type;
        typedef t_rank                               rank_1_type;
        typedef t_select1                            select_1_type;
        typedef t_select0                            select_0_type;
        typedef n_bitvector                          nbit_vector_type;
        typedef n_rank                               nrank_1_type;
        typedef n_select1                            nselect_1_type;
        typedef n_select0                            nselect_0_type;
        typedef t_succ_tree                          succ_tree;

protected:
wt_int<> wt1,wtd1,wtnd1,wtns1;
wt_int<> wtd2,wtnd2,wtns2;
wt_int<> wt2,wt2maped,wt2mapedS;
int_vector<> ids, ids2, ids3, ids4   , dids1,dids2,dids3;
vector<int> original,original2,original3,original4;
vector <bool> inicial,inicial2,inicial3,inicial4;
vector <int> B,B2,B3,B4,H,H2,H3,H4,HM;
inv_perm_support<> perm, perm2, perm3, perm4;
int mitotal = 0;
bp_support_sada<256,32,rank_support_v<>,select_support_mcl<>> bps,bpsd,bpsnd,bpsns;
//bp2 bitvector marca 1 los 2, bopen marca 1 los 1
bit_vector Baux2, bp, bopen,bp2,hoja, bopend, bopennd, bopenns,bp2d,bp2nd, bp2ns,bp22,bp3,bp4,mapbarato,bp2H;
sd_vector<> bp2maped;
sd_vector<>::select_1_type bp2maped_select;
bit_vector::rank_1_type b_rank, b2_rank,bopen_rank,bopend_rank, bopennd_rank, bopenns_rank,bp2d_rank,bp2nd_rank, bp2ns_rank,b2H_rank;
bit_vector::select_1_type b_select, b2_select,bopen_select,bopend_select, bopennd_select, bopenns_select,bp2d_select,bp2nd_select, bp2ns_select;//,bp2maped_select;
vector <int> wtbarato, permbarato;
int limitenodos,maxwt2,minwt2;

        void copy(const gbp& p) {
            original          = p.original;
            original2         = p.original2;
            original3 		  = p.original3;
            original4         = p.original4;
            inicial          = p.inicial;
            inicial2         = p.inicial2;
            inicial3        = p.inicial3;
            inicial4       = p.inicial4;
            ids       = p.ids;
            ids2       = p.ids2;
            ids3       = p.ids3;
            ids4  = p.ids4;
            wtbarato  = p.wtbarato;
            dids1  = p.dids1;
            dids2     = p.dids2;
            dids3     = p.dids3;
            B       = p.B;
            B2  = p.B2;
            B3     = p.B3;
            B4     = p.B4;
            H     = p.H;
            H2  = p.H2;
            H3     = p.H3;
            H4  = p.H4;
            perm     = p.perm;
            perm2     = p.perm2;
            perm3     = p.perm3;
            perm4     = p.perm4;
            limitenodos		= p.limitenodos;
        	maxwt2			= p.maxwt2;
        	minwt2			= p.minwt2;
            HM  = p.HM;        	
            Baux2 = p.Baux2;
            bp = p.bp; 
            bopen = p.bopen;
            bp2 = p.bp2;
            hoja = p.hoja; 
            bopend = p.bopend; 
            bopennd = p.bopennd; 
            bopenns = p.bopenns;
            bp2d = p.bp2d;
            bp2nd = p.bp2nd; 
            bp2ns = p.bp2ns;
            bp22 = p.bp22;
            bp3 = p.bp3;
            bp4 = p.bp4;
            mapbarato = p.mapbarato;
            bp2H = p.bp2H;
            wtbarato = p.wtbarato; 
            permbarato = p.permbarato;
            bp2maped = p.bp2maped;
        }


    public:

        //! Default constructor
        gbp() {};
//sep default = 1;
     	gbp(string file1, string file2, string file3, string file4) {
     int granularidades, aux,mitotal = 0,relaciones,a,b,c,d;
FILE *fp = fopen(file1.c_str(), "r");
    if (!fp) {
    fprintf(stderr, "Error opening file \"%s\".\n", file1.c_str());
    exit(EXIT_FAILURE);
    }
   fscanf(fp,"%d",&granularidades);
    vector <int> acumulado,totalnivel(granularidades,0);

    for(int i = 0; i < granularidades; i++){
        fscanf(fp,"%d",&aux);
        acumulado.push_back(mitotal);
        mitotal+= aux;
        totalnivel[i] = aux;
    }
    limitenodos = mitotal;
    fscanf(fp,"%d",&relaciones);
    vector <vector<int> > adjlist(mitotal, vector<int>());
    for(int i =  0; i < relaciones; i++){
        fscanf(fp,"%d %d %d %d",&a,&b,&c,&d);
        adjlist[c+acumulado[d]].push_back(a+acumulado[b]);
   }

    ids.resize(mitotal);
    util::set_to_value(ids,0);
    int posactual = 0;
    vector <int> visitados(mitotal, false);
    vector <int> cerrar(mitotal, false);
    for(int i = 0; i < totalnivel[granularidades-1]; i++){

        int actual = i + acumulado[granularidades-1];
        stack <int> q;
        q.push(actual);
        while(!q.empty()){
            int frente = q.top();
            if(cerrar[frente]){
                q.pop();
                B.push_back(0);
                cerrar[frente] = false;
                continue;
            }
            else {
                cerrar[frente] = true;
                if(!visitados[frente]){
                    B.push_back(1);
                    ids[frente] = posactual;
                    posactual++;
                } 
                else{
                    B.push_back(2);
                    H.push_back(ids[frente]);
                    visitados[frente] = true;
                    continue;
                } 
            } 
            visitados[frente] = true;
            for(int j = 0; j < adjlist[frente].size();j++){
                    if(!cerrar[adjlist[frente][j]])q.push(adjlist[frente][j]);
            }

        }
    }
int contadorabiertos = 0;
    bit_vector auxA(B.size(),0),auxB(B.size(),0), auxC(B.size(),0),auxD(B.size(),0),auxE(B.size(),0),auxF(B.size(),0),auxJ(B.size()+1,0),auxK(B.size(),0),auxI(B.size()/2,0);
for(int i = 0 ; i < B.size(); i++){
    if(B[i] > 0)auxC[i] = 1;
    if(B[i] == 1)auxA[i] = 1;
    else if(B[i] == 2)auxB[i] = 1;
    if(B[i] == 2)auxD[i] = 1;
    if(B[i] == 1)auxE[i] = 1;
if(B[i] == 2)auxI[contadorabiertos] = 1;
    if(B[i] >= 1)contadorabiertos++;
}
int totalhojas = 0;
for(int i = 0; i < B.size()-1; i++)if(B[i] >= 1 && B[i+1] == 0){
    auxF[i] = 1;
    totalhojas++;
}

bp_support_sada<256,32,rank_support_v<>,select_support_mcl<>> bpsaux(&auxC);
bps.swap(bpsaux);
bp.swap(auxC);
bopen.swap(auxE);
bps.set_vector(&bp);
bp2.swap(auxD);
bp2H.swap(auxI);
util::init_support(b_rank, &auxD);
util::init_support(b2_rank, &bp2);
util::init_support(b2H_rank, &bp2H);
util::init_support(b_select, &auxD);
util::init_support(b2_select, &bp2);
util::init_support(bopen_select, &bopen);
util::init_support(bopen_rank, &bopen);



for(int i = 0; i < H.size();i++)HM.push_back(bopen_select(H[i]+1));


int_vector<> vB(B.size());
int_vector <> vH(H.size());
int_vector <> vHM(H.size());
int_vector <> vHM2(H.size());
map <int,int> mapcont;
map <int,int> mapcont2;
    for(int i = 0 ; i < B.size(); i++)vB[i] = B[i];
    for(int i = 0 ; i < H.size(); i++)vH[i] = H[i];
for(int i = 0 ; i < HM.size(); i++)vHM2[i] = HM[i];
	for(int i = 0; i < H.size();i++)mapcont[H[i]] = 1;
int contadormapeos = 0;
for(int i = 0; i < B.size()+1;i++){
	if(mapcont[i] == 1){
	auxJ[i] = 1;
	mapcont2[i] = contadormapeos;
	contadormapeos++;
	}
}
for(int i = 0; i < H.size(); i++)vHM[i] = mapcont2[H[i]];

sd_vector <> auxJJ(auxJ);
bp2maped.swap(auxJJ);


util::init_support(bp2maped_select, &bp2maped); 
construct_im(wt1,vB);
construct_im(wt2,vH);
construct_im(wt2maped,vHM);
construct_im(wt2mapedS,vHM2);
inv_perm_support<> auxperm1(&ids);
perm = auxperm1;

bit_vector auxG(wt2.size(),0);
int posinb2 = 0;
totalhojas = 0;
for(int  i = 0; i < B.size();i++){
    if(auxF[i] == 1 && bp2[i] == 1){
        auxG[posinb2] = 1;
        totalhojas++;
    }
    if(bp2[i]==1)posinb2++;
}
hoja.swap(auxG);


for(int i = 0; i < wt2mapedS.size();i++)if(hoja[i] == 1)auxK[wt2mapedS[i]] = 1;
mapbarato.swap(auxK);
maxwt2 = 0;
for(int i = 0; i < wt2.size();i++)if(wt2[i] > maxwt2)maxwt2 = wt2[i];
minwt2 = maxwt2;
for(int i = 0; i < wt2.size();i++)if(wt2[i] < minwt2)minwt2 = wt2[i];
fclose(fp);

    FILE *fp2 = fopen(file2.c_str(), "r");
    if (!fp2) {
    fprintf(stderr, "Error opening file \"%s\".\n", file2.c_str());
    exit(EXIT_FAILURE);
    }
      fscanf(fp2,"%d",&relaciones);
    vector <vector<int> > adjlist2(mitotal+1, vector<int>());
    map <pair<int,int>,int> mapita2;
    for(int i =  0; i < relaciones; i++){
        fscanf(fp2,"%d %d %d %d",&a,&b,&c,&d);
        if(mapita2[make_pair(a+acumulado[b],c+acumulado[d])]== 1) continue;
        mapita2[make_pair(a+acumulado[b],c+acumulado[d])]++;
        mapita2[make_pair(c+acumulado[d],a+acumulado[b])]++;
        adjlist2[a+acumulado[b]].push_back(c+acumulado[d]);
        adjlist2[c+acumulado[d]].push_back(a+acumulado[b]);
    }
    ids2.resize(mitotal+1);
    util::set_to_value(ids2,0);
    posactual = 1;
    vector <bool> visitados2(mitotal+1, false),unidos(mitotal+1,false);

    for(int i = mitotal-1; i >= 0; i--){
        if(!unidos[i]){
            unidos[i] = true;
            adjlist2[mitotal].push_back(i);
            queue<int> q;
            q.push(i);
            while(!q.empty()){
                int frente = q.front();
                q.pop();
                for(int j = 0;j < adjlist2[frente].size();j++){
                    if(!unidos[adjlist2[frente][j]]){
                        unidos[adjlist2[frente][j]] = true;
                        q.push(adjlist2[frente][j]);
                    }
                }
            }
        }
    }

        stack <int> q,st;
    visitados2[mitotal] = true;
    vector <int> cerrar2(mitotal+1, false);
    for(int i = 0; i < adjlist2[mitotal].size(); i++){
                int actual = adjlist2[mitotal][i];
                st.push(actual);
    }


        while(!st.empty()){
            int frente = st.top();
            if(cerrar2[frente]){
                st.pop();
                B2.push_back(0);
                cerrar2[frente] = false;
                continue;
            }
            else {
                cerrar2[frente] = true;
                if(!visitados2[frente]){
                    B2.push_back(1);
                    ids2[frente] = posactual;
                    posactual++;
                } 
                else{
                    B2.push_back(2);
                    H2.push_back(ids2[frente]);
                    visitados2[frente] = true;
                    continue;
                } 
            } 
            visitados2[frente] = true;
            for(int j = 0; j < adjlist2[frente].size();j++){
                    if(!cerrar2[adjlist2[frente][j]])st.push(adjlist2[frente][j]);
            }

        }
    fclose(fp2);

    FILE *fp3 = fopen(file3.c_str(), "r");
    if (!fp3) {
    fprintf(stderr, "Error opening file \"%s\".\n", file3.c_str());
    exit(EXIT_FAILURE);
    }
    fscanf(fp3,"%d",&relaciones);
            map <pair<int,int>,int> comprepe;
    vector <vector<int> > adjlist3(mitotal+1, vector<int>());
    for(int i =  0; i < relaciones; i++){
        fscanf(fp3,"%d %d %d %d",&a,&b,&c,&d);
             if(comprepe[make_pair(a+acumulado[b],c+acumulado[d])]== 1) continue;
                                        comprepe[make_pair(a+acumulado[b],c+acumulado[d])]++;
                                comprepe[make_pair(c+acumulado[d],a+acumulado[b])]++;
        adjlist3[a+acumulado[b]].push_back(c+acumulado[d]);
        adjlist3[c+acumulado[d]].push_back(a+acumulado[b]);
    }

    ids3.resize(mitotal+1);
    util::set_to_value(ids3,0);
    posactual = 1;
    visitados2.clear();
    unidos.clear();
    visitados2.resize(mitotal+1,false);
    unidos.resize(mitotal+1,false);
    stack <int> q2;
    for(int i = mitotal-1; i >= 0; i--){
        if(!unidos[i]){
            unidos[i] = true;
            adjlist3[mitotal].push_back(i);
            queue<int> q;
            q.push(i);
            while(!q.empty()){
                int frente = q.front();
                q.pop();
                for(int j = 0;j < adjlist3[frente].size();j++){
                    if(!unidos[adjlist3[frente][j]]){
                        unidos[adjlist3[frente][j]] = true;
                        q.push(adjlist3[frente][j]);
                    }
                }
            }
        }
    }

    visitados2[mitotal] = true;
    for(int i = 0; i < adjlist3[mitotal].size(); i++){
                int actual = adjlist3[mitotal][i];
               q2.push(actual);
    }

    vector <int> cerrar3(mitotal, false);
        while(!q2.empty()){
            int frente = q2.top();
            if(cerrar3[frente]){
                q2.pop();
                B3.push_back(0);
                cerrar3[frente] = false;
                continue;
            }
            else {
                cerrar3[frente] = true;
                if(!visitados2[frente]){
                    B3.push_back(1);
                    ids3[frente] = posactual;
                    posactual++;
                } 
                else{
                    B3.push_back(2);
                    H3.push_back(ids3[frente]);
                    visitados2[frente] = true;
                    continue;
                } 
            } 
            visitados2[frente] = true;
            for(int j = 0; j < adjlist3[frente].size();j++){
                    if(!cerrar3[adjlist3[frente][j]])q2.push(adjlist3[frente][j]);
            }

        }
    fclose(fp3);


 FILE *fp4 = fopen(file4.c_str(), "r");
    if (!fp4) {
    fprintf(stderr, "Error opening file \"%s\".\n",file4.c_str());
    exit(EXIT_FAILURE);
    }
    fscanf(fp4,"%d",&relaciones);
    vector <vector<int> > adjlist4(mitotal+1, vector<int>());
    map <pair<int,int>,int> mapita4;
    for(int i =  0; i < relaciones; i++){
        fscanf(fp4,"%d %d %d %d",&a,&b,&c,&d);
        if(mapita4[make_pair(a+acumulado[b],c+acumulado[d])]== 1) continue;
        mapita4[make_pair(a+acumulado[b],c+acumulado[d])]++;
        mapita4[make_pair(c+acumulado[d],a+acumulado[b])]++;
        adjlist4[a+acumulado[b]].push_back(c+acumulado[d]);
        adjlist4[c+acumulado[d]].push_back(a+acumulado[b]);
    }

    ids4.resize(mitotal+1);
    util::set_to_value(ids4,0);
    posactual = 1;
    visitados2.clear();
    unidos.clear();
    visitados2.resize(mitotal+1,false);
    unidos.resize(mitotal+1,false);

    for(int i = mitotal-1; i >= 0; i--){
        if(!unidos[i]){
            unidos[i] = true;
            adjlist4[mitotal].push_back(i);
            queue<int> q;
            q.push(i);
            while(!q.empty()){
                int frente = q.front();
                q.pop();
                for(int j = 0;j < adjlist4[frente].size();j++){
                    if(!unidos[adjlist4[frente][j]]){
                        unidos[adjlist4[frente][j]] = true;
                        q.push(adjlist4[frente][j]);
                    }
                }
            }
        }
    }
    visitados2[mitotal] = true;
    stack <int> sq3;
    for(int i = 0; i < adjlist4[mitotal].size(); i++){
                int actual = adjlist4[mitotal][i];
                sq3.push(actual);
    }
//    cout << "tam de queue: " << adjlist4[mitotal].size() << endl;
    vector <int> cerrar4(mitotal+1, false);
        while(!sq3.empty()){
            int frente = sq3.top();
            if(cerrar4[frente]){
                sq3.pop();
                B4.push_back(0);
                cerrar4[frente] = false;
                continue;
            }
            else {
                cerrar4[frente] = true;
                if(!visitados2[frente]){
                    B4.push_back(1);
                    ids4[frente] = posactual;
                    posactual++;
                } 
                else{
                    B4.push_back(2);
                    H4.push_back(ids4[frente]);
                    visitados2[frente] = true;
                    continue;
                } 
            } 
            visitados2[frente] = true;
            for(int j = 0; j < adjlist4[frente].size();j++){
                    if(!cerrar4[adjlist4[frente][j]])sq3.push(adjlist4[frente][j]);
            }

        }
    fclose(fp4);


int auxcant1 = 0;
for(int i = 0; i < B3.size(); i++)if(B3[i]>=1)auxcant1++;


    int_vector<> vB2(B2.size()),vB3(B3.size()) ,vB4(B4.size()); 
    int_vector <> vH2(H2.size()),vH3(H3.size()), vH4(H4.size()); 
    sd_vector <> prueba1(auxA);
    sd_vector <> prueba2(auxB);
    rrr_vector <63> prueba3(auxA);
    rrr_vector <63> prueba4(auxB);


    for(int i = 0 ; i < B2.size(); i++)vB2[i] = B2[i];
    for(int i = 0 ; i < H2.size(); i++)vH2[i] = H2[i];

    for(int i = 0 ; i < B3.size(); i++)vB3[i] = B3[i];
    for(int i = 0 ; i < H3.size(); i++)vH3[i] = H3[i];

    for(int i = 0 ; i < B4.size(); i++)vB4[i] = B4[i];
    for(int i = 0 ; i < H4.size(); i++)vH4[i] = H4[i];


construct_im(wtd1,vB2);
construct_im(wtd2,vH2);
//not subsumption
construct_im(wtns1,vB3);
construct_im(wtns2,vH3);
//not disjoint
construct_im(wtnd1,vB4);
construct_im(wtnd2,vH4);


//INIT RESTO DE GLOUDBP
bit_vector auxA2(B2.size(),0),auxB2(B2.size(),0), auxC2(B2.size(),0);
for(int i = 0 ; i < B2.size(); i++){
    if(B2[i] > 0)auxA2[i] = 1;
    if(B2[i] == 1)auxB2[i] = 1;
    else if(B2[i] == 2)auxC2[i] = 1;
}
bp_support_sada<256,32,rank_support_v<>,select_support_mcl<>> bpsaux2(&auxA2);
bopend.swap(auxB2);
bp2d.swap(auxC2);
bpsd.swap(bpsaux2);
util::init_support(bp2d_rank, &bp2d);
util::init_support(bp2d_select, &bp2d);
util::init_support(bopend_select, &bopend);
util::init_support(bopend_rank, &bopend);


bit_vector auxA3(B3.size(),0),auxA33(B3.size(),0),auxB3(B3.size(),0), auxC3(B3.size(),0);
for(int i = 0 ; i < B3.size(); i++){
    if(B3[i] > 0)auxA3[i] = 1;
    if(B3[i] > 0)auxA33[i] = 1;
    if(B3[i] == 1)auxB3[i] = 1;
    else if(B3[i] == 2)auxC3[i] = 1;
}
bp_support_sada<256,32,rank_support_v<>,select_support_mcl<>> bpsaux3(&auxA3);
bopenns.swap(auxB3);
bp2ns.swap(auxC3);
bpsns.swap(bpsaux3);
util::init_support(bp2ns_rank, &bp2ns);
util::init_support(bp2ns_select, &bp2ns);
util::init_support(bopenns_select, &bopenns);
util::init_support(bopenns_rank, &bopenns);



bit_vector auxA4(B4.size(),0),auxB4(B4.size(),0), auxC4(B4.size(),0);
for(int i = 0 ; i < B4.size(); i++){
    if(B4[i] > 0)auxA4[i] = 1;
    if(B4[i] == 1)auxB4[i] = 1;
    else if(B4[i] == 2)auxC4[i] = 1;
}
bp_support_sada<256,32,rank_support_v<>,select_support_mcl<>> bpsaux4(&auxA4);
bopennd.swap(auxB4);
bp2nd.swap(auxC4);
bpsnd.swap(bpsaux4);
util::init_support(bp2nd_rank, &bp2nd);
util::init_support(bp2nd_select, &bp2nd);
util::init_support(bopennd_select, &bopennd);
util::init_support(bopennd_rank, &bopennd);





inv_perm_support<> auxperm2(&ids2);
inv_perm_support<> auxperm3(&ids3);
inv_perm_support<> auxperm4(&ids4);


perm2 = auxperm2;
perm3 = auxperm3;
perm4 = auxperm4;

bp22.swap(auxA2);
bp3.swap(auxA33);
bp4.swap(auxA4);
bpsd.set_vector(&bp22);
bpsns.set_vector(&bp3);
bpsnd.set_vector(&bp4);

util::bit_compress(ids);
util::bit_compress(ids2);
util::bit_compress(ids3);
util::bit_compress(ids4);

for(int i = 0; i < wt2.size();i++)wtbarato.push_back(wt2[i]);

map<int,int>repetidos;
int contarepetidos = 0;
for(int i = 0; i < wtbarato.size();i++){
if(repetidos[wtbarato[i]]== 0){
    repetidos[wtbarato[i]]++;
    contarepetidos++;
}
}


for(int i = 0 ; i < ids.size(); i++) permbarato.push_back(perm[i]);


    dids1.resize(mitotal);
    util::set_to_value(dids1,0);
    dids2.resize(mitotal);
    util::set_to_value(dids2,0);
    dids3.resize(mitotal);
    util::set_to_value(dids3,0);
    for(int i = 0; i < mitotal; i++){
        dids1[i] = mapOtoB(mapBtoO(i,0),1);
    }
    for(int i = 0; i < mitotal; i++){
        dids2[i] = mapOtoB(mapBtoO(i,0),2);
    }
    for(int i = 0; i < mitotal; i++){
        dids3[i] = mapOtoB(mapBtoO(i,0),3);
    }

util::bit_compress(dids1);
util::bit_compress(dids2);
util::bit_compress(dids3);

 }

        //! Copy constructor
        gbp(const gbp& g) {
            copy(g);
        }

        //! Copy constructor
        gbp(gbp&& g) {
            *this = std::move(g);
        }

        //! Assignment operator
        gbp& operator=(const gbp g) {
            if (this != &g) {
                copy(g);
            }
            return *this;
        }

        //! Assignment move operator
        gbp& operator=(gbp&& g) {
            if (this != &g) {
            original          = std::move(g.original);
            original2         = std::move(g.original2);
            original3 		  = std::move(g.original3);
            original4         = std::move(g.original4);
            inicial          = std::move(g.inicial);
            inicial2         = std::move(g.inicial2);
            inicial3        = std::move(g.inicial3);
            inicial4       = std::move(g.inicial4);
            ids       = std::move(g.ids);
            ids2       = std::move(g.ids2);
            ids3       = std::move(g.ids3);
            ids4  = std::move(g.ids4);
            wtbarato  = std::move(g.wtbarato);
            dids1  = std::move(g.dids1);
            dids2     = std::move(g.dids2);
            dids3     = std::move(g.dids3);
            B       = std::move(g.B);
            B2  = std::move(g.B2);
            B3     = std::move(g.B3);
            B4     = std::move(g.B4);
            H     = std::move(g.H);
            H2  = std::move(g.H2);
            H3     = std::move(g.H3);
            H4  = std::move(g.H4);
            perm     = std::move(g.perm);
            perm2     = std::move(g.perm2);
            perm3     = std::move(g.perm3);
            perm4     = std::move(g.perm4);

            HM  = std::move(g.HM);             
            Baux2 = std::move(g.Baux2);
            bp = std::move(g.bp); 
            bopen = std::move(g.bopen);
            bp2 = std::move(g.bp2);
            hoja = std::move(g.hoja); 
            bopend = std::move(g.bopend); 
            bopennd = std::move(g.bopennd); 
            bopenns = std::move(g.bopenns);
            bp2d = std::move(g.bp2d);
            bp2nd = std::move(g.bp2nd); 
            bp2ns = std::move(g.bp2ns);
            bp22 = std::move(g.bp22);
            bp3 = std::move(g.bp3);
            bp4 = std::move(g.bp4);
            mapbarato = std::move(g.mapbarato);
            bp2H = std::move(g.bp2H);
            wtbarato = std::move(g.wtbarato); 
            permbarato = std::move(g.permbarato);
            bp2maped = std::move(g.bp2maped);
            }
            return *this;
        }

        //! Swap operator
        void swap(gbp& g) {/*
            if (this != &g) {
                std::swap(m_vertices, g.m_vertices);
                std::swap(m_levels, g.m_levels);
                std::swap(m_edges,  g.m_edges);
		
		m_A.swap(g.m_A);
                util::swap_support(m_A_rank, g.m_A_rank, &m_A, &(g.m_A));
                util::swap_support(m_A_select1, g.m_A_select1, &m_A, &(g.m_A));
                util::swap_support(m_A_select0, g.m_A_select0, &m_A, &(g.m_A));
		B.swap(g.B);
                util::swap_support(B_rank, g.B_rank, &B, &(g.B));
                util::swap_support(B_select1, g.B_select1, &B, &(g.B));

		BP.swap(g.BP);
                util::swap_support(BP_rank, g.BP_rank, &BP, &(g.BP));
                util::swap_support(BP_select1, g.BP_select1, &BP, &(g.BP));  
        BF.swap(g.BF);
        		util::swap_support(BF_rank, g.BF_rank, &BF, &(g.BF));
        		util::swap_support(BF_select1, g.BF_select1, &BF, &(g.BF));
        		util::swap_support(BF_select0, g.BF_select0, &BF, &(g.BF));
        		util::swap_support(BF_st, g.BF_st, &BF, &(g.BF));

		B_F.swap(g.B_F);
                util::swap_support(B_F_rank, g.B_F_rank, &B_F, &(g.B_F));
                util::swap_support(B_F_select1, g.B_F_select1, &B_F, &(g.B_F));

		m_B.swap(g.m_B);
		        util::swap_support(m_B_rank, g.m_B_rank, &m_B, &(g.m_B));
                util::swap_support(m_B_select1, g.m_B_select1, &m_B, &(g.m_B));
                util::swap_support(m_B_st, g.m_B_st, &m_B, &(g.m_B));
                util::swap_support(m_B_select0, g.m_B_select0, &m_B, &(g.m_B));

		m_B_star.swap(g.m_B_star);
                util::swap_support(m_B_star_st, g.m_B_star_st, &m_B_star,
				   &(g.m_B_star));
            }
            */
        }



        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
	  structure_tree_node* child = structure_tree::add_child(v, name,
								 util::class_name(*this));

	  size_type written_bytes = 0;



	  written_bytes += wtns1.serialize(out, child, "wtns1");
	  written_bytes += wt1.serialize(out, child, "wt1");
	  written_bytes += wtnd1.serialize(out, child, "wtnd1");
	  written_bytes += wtd1.serialize(out, child, "wtd1");
	  written_bytes += bps.serialize(out, child, "bps");
	  written_bytes += bpsd.serialize(out, child, "bpsd");
	  written_bytes += bpsnd.serialize(out, child, "bpsnd");
	  written_bytes += bpsns.serialize(out, child, "bpsns");

	  written_bytes += bopend.serialize(out, child, "bopend");
	  written_bytes += bopennd.serialize(out, child, "bopennd");
	  written_bytes += bopenns.serialize(out, child, "bopenns");
	  written_bytes += bopend_rank.serialize(out, child, "bopend_rank");
	  written_bytes += bopennd_rank.serialize(out, child, "bopennd_rank");
	  written_bytes += bopenns_rank.serialize(out, child, "bopenns_rank");
	  written_bytes += bopend_select.serialize(out, child, "bopend_select");
	  written_bytes += bopennd_select.serialize(out, child, "bopennd_select");
	  written_bytes += bopenns_select.serialize(out, child, "bopenns_select");

	  written_bytes += wtns2.serialize(out, child, "wtns2");
	  written_bytes += wtnd2.serialize(out, child, "wtnd2");
	  written_bytes += wt2mapedS.serialize(out, child, "wt2mapedS");
	  written_bytes += wtd2.serialize(out, child, "wtd2");

	  written_bytes += bp2d.serialize(out, child, "bp2d");
	  written_bytes += bp2nd.serialize(out, child, "bp2nd");
	  written_bytes += bp2ns.serialize(out, child, "bp2ns");
	  written_bytes += bp2d_rank.serialize(out, child, "bp2d_rank");
	  written_bytes += bp2nd_rank.serialize(out, child, "bp2nd_rank");
	  written_bytes += bp2ns_rank.serialize(out, child, "bp2ns_rank");
	  written_bytes += bp2d_select.serialize(out, child, "bp2d_select");
	  written_bytes += bp2nd_select.serialize(out, child, "bp2nd_select");
	  written_bytes += bp2ns_select.serialize(out, child, "bp2ns_select");

	  written_bytes += bp.serialize(out, child, "bp");
	  written_bytes += bopen.serialize(out, child, "bopen");
	  written_bytes += bp2H.serialize(out, child, "bp2H");
	  written_bytes += bps.serialize(out, child, "bps");
	  written_bytes += b2H_rank.serialize(out, child, "b2H_rank");
	  written_bytes += bopen_rank.serialize(out, child, "bopen_rank");
	  written_bytes += bopen_select.serialize(out, child, "bopen_select");


	  written_bytes += perm.serialize(out, child, "perm");
	  written_bytes += ids.serialize(out, child, "ids");
	  written_bytes += dids1.serialize(out, child, "dids1");
	  written_bytes += dids2.serialize(out, child, "dids2");
	  written_bytes += dids3.serialize(out, child, "dids3");

	  structure_tree::add_size(child, written_bytes);
	  return written_bytes;
        }

        //! Loads the data structure from the given istream.




        void load(std::istream& in) {
	    wt1.load(in);
	    wtnd1.load(in);
	    wtd1.load(in);

	    wtns2.load(in);
	    wtnd2.load(in);
	    wt2.load(in);
	    wtd2.load(in);

	    perm.load(in);
	    ids.load(in);
	    dids1.load(in);
	    dids2.load(in);
	    dids3.load(in);
          
            Baux2.load(in);
            bp.load(in); 
            bopen.load(in);
            bp2.load(in);
            hoja.load(in); 
            bopend.load(in); 
            bopennd.load(in); 
            bopenns.load(in);
            bp2d.load(in);
            bp2nd.load(in); 
            bp2ns.load(in);
            bp22.load(in);
            bp3.load(in);
            bp4.load(in);
            mapbarato.load(in);
            bp2H.load(in);
            bp2maped.load(in);
        }


//RECORDAR QUE PARA TIPO = 2,3,4. SE PARTE DESDE 1 DADO QUE 0 ES NODO FICTICIO CREADO

int mapbp2(int id, int tipo){
    if(tipo == 0) return bp2[id];
    if(tipo == 1) return bp2d[id];
    if(tipo == 2) return bp2ns[id];
    if(tipo == 3) return bp2nd[id];
    return -1;
}


int mapBtoO(int id, int tipo){
    if(tipo == 0) return perm[id];
    if(tipo == 1) return perm2[id];
    if(tipo == 2) return perm3[id];
    if(tipo == 3) return perm4[id];
    return -1;
}

/*
Dado un identificador en O, retornar el identificador mapeado a indices en B
tipo refiere a cual de los gloud, sinedo 1= contencion, 2 = disjoint, 3 = nocontencion, 4 = notdisjoint
*/
int mapOtoB(int id, int tipo){
    if(tipo == 0) return ids[id];
    if(tipo == 1) return ids2[id];
    if(tipo == 2) return ids3[id];
    if(tipo == 3) return ids4[id];
    return -1;
}


vector <int> padres(int id){
/*
Funcion para obtener todos los padres dado un nodo dada un id en B
Si quiero unicamente los padres directo simplemente quito el while (que itere unicamente 1 vez), dado que
un nodo en un arbol tiene unicamente un padre sera un retorno directo
*/
vector <int> papas;
if(id == 0)return papas;
map<int,int> comprobar;
queue<int> q;
//Obtengo su posicion en el wt que representa la secuencia B
// con esta linea mapeo de posicion en B a identificador en B
// para mapear de identificador a posicion se realiza un select(identificador,1)
comprobar[id] = 1;
q.push(id);
    while(!q.empty()){
    id = q.front();
    q.pop();
    int nuevoB = wt1.rank(wt1.select(id,1),0);
    if(comprobar[nuevoB ] == 0 ){
        papas.push_back(nuevoB );
        if(!inicial[nuevoB])q.push(nuevoB);
        comprobar[nuevoB]++;
    }
    //Limite se calcula asi: busco si para el nodo actualmente analizado existe alguna aparicion en H, de ser asi obtengo todas su apariciones
    int j = 1, limite = wt2.rank(wt2.size(),id);
    while(j <= limite){
        nuevoB = wt1.rank(wt1.select(wt2.select(j,id)+1,2),0);
      //  cout << "entre con id " << id <<endl;
        if(comprobar[nuevoB] == 0 ){
            papas.push_back(nuevoB);
            if(!inicial[nuevoB])q.push(nuevoB);
            comprobar[nuevoB]++;
        }
        j++;
    }
    }
    return papas;

}


//ID es identificador de un nodo en B, no su posicion
void hijo(int id){
 vector <int> sons;
map<int,int> comprobar;
queue<int> q;   
//sumo 1 porque al ser un select e identificador partir de 0 debo sumar un 1.
//Obtengo la posicion donde terminan los hijos del nodo id
id = wt1.select(id+1,0);
comprobar[id] = 1;
q.push(id);
for(int i = id-1; i >=0; i--){
    if(wt1[i] == 0)break;
    if(wt1[i] == 1)sons.push_back(wt1.rank(i+1,1));
    else sons.push_back(wt2[wt1.rank(i+1,2)-1]);
}
}


int mapeogetrankt(int id, int tipo){
    if(tipo == 0) return bps.rank(id);
    if(tipo == 1) return bpsd.rank(id);
    if(tipo == 2) return bpsns.rank(id);
    if(tipo == 3) return bpsnd.rank(id);
    return -1;
}

int mapeogetselect(int id, int tipo){
    if(tipo == 0) return bps.rank(id);
    if(tipo == 1) return bpsd.rank(id);
    if(tipo == 2) return bpsns.rank(id);
    if(tipo == 3) return bpsnd.rank(id);
    return -1;
}

//ID es identificador de un nodo en B, no su posicion
int mapeo1select(int id, int tipo, int x){
    if(tipo == 0) return wt1.select(id,x);
    if(tipo == 1) return wtd1.select(id,x);
    if(tipo == 2) return wtns1.select(id,x);
    if(tipo == 3) return wtnd1.select(id,x);
    return -1;
}

int mapeogetclose(int id, int tipo){
    if(tipo == 0) return bps.find_close(id);
    if(tipo == 1) return bpsd.find_close(id);
    if(tipo == 2) return bpsns.find_close(id);
    if(tipo == 3) return bpsnd.find_close(id);
    return -1;
}


int mapeogetenclose(int id, int tipo){
    if(tipo == 0) return bps.enclose(id);
    if(tipo == 1) return bpsd.enclose(id);
    if(tipo == 2) return bpsns.enclose(id);
    if(tipo == 3) return bpsnd.enclose(id);
    return -1;
}


int mapeo1rankbopen(int id, int tipo){
    if(tipo == 0) return bopen_rank(id+1);
    if(tipo == 1) return bopend_rank(id+1);
    if(tipo == 2) return bopenns_rank(id+1);
    if(tipo == 3) return bopennd_rank(id+1);
    return -1;
}


int mapeorankbp2(int id, int tipo){
    if(tipo == 0) return b2_rank(id+1);
    if(tipo == 1) return bp2d_rank(id+1);
    if(tipo == 2) return bp2ns_rank(id+1);
    if(tipo == 3) return bp2nd_rank(id+1);
    return -1;
}


int mapeo1selectbopen(int id, int tipo){
    if(tipo == 0) return bopen_select(id+1);
    if(tipo == 1) return bopend_select(id+1);
    if(tipo == 2) return bopenns_select(id+1);
    if(tipo == 3) return bopennd_select(id+1);
    return -1;
}


int mapeoselectbp2(int id, int tipo){
    if(tipo == 0) return b2_select(id+1);
    if(tipo == 1) return bp2d_select(id+1);
    if(tipo == 2) return bp2ns_select(id+1);
    if(tipo == 3) return bp2nd_select(id+1);
    return -1;
}

int mapeo1rank(int id, int tipo, int x){
    if(tipo == 0) return wt1.rank(id,x);
    if(tipo == 1) return wtd1.rank(id,x);
    if(tipo == 2) return wtns1.rank(id,x);
    if(tipo == 3) return wtnd1.rank(id,x);
    return -1;
}

//ID es identificador de un nodo en B, no su posicion
int mapeo2select(int id, int tipo, int x){
    if(tipo == 0) return wt2.select(id,x);
    if(tipo == 1) return wtd2.select(id,x);
    if(tipo == 2) return wtns2.select(id,x);
    if(tipo == 3) return wtnd2.select(id,x);
    return -1;
}


int mapeo2rank(int id, int tipo, int x){
    if(tipo == 0) return wt2.rank(id,x);
    if(tipo == 1) return wtd2.rank(id,x);
    if(tipo == 2) return wtns2.rank(id,x);
    if(tipo == 3) return wtnd2.rank(id,x);
    return -1;
}


//POSICION A EVALUAR EN B, NO ID
int mapeodirecto1(int id, int tipo){
    if(tipo == 0) return wt1[id];
    if(tipo == 1) return wtd1[id];
    if(tipo == 2) return wtns1[id];
    if(tipo == 3) return wtnd1[id];
    return -1;
}

//POSICION A EVALUAR EN H, NO ID
int mapeodirecto2(int id, int tipo){
    if(tipo == 0) return wt2[id];
    if(tipo == 1) return wtd2[id];
    if(tipo == 2) return wtns2[id];
    if(tipo == 3) return wtnd2[id];
    return -1;
}


 int tam3(int tipo){
    if(tipo == 0) return bps.size();
    if(tipo == 1) return bpsd.size();
    if(tipo == 2) return bpsns.size();
    if(tipo == 3) return bpsnd.size();
 }


 int tam2(int tipo){
    if(tipo == 0) return wt2.size();
    if(tipo == 1) return wtd2.size();
    if(tipo == 2) return wtns2.size();
    if(tipo == 3) return wtnd2.size();
 }

 int tam1(int tipo){
    if(tipo == 0) return wt1.size();
    if(tipo == 1) return wtd1.size();
    if(tipo == 2) return wtns1.size();
    if(tipo == 3) return wtnd1.size();
 }


//ID es identificador de un nodo en B, no su posicion
vector <int> vecinos(int id,int tipo){
 vector <int> sons;
map<int,int> comprobar,agregados;
queue<int> q;   
int idoriginal = id;

//sumo 1 porque al ser un select e identificador partir de 0 debo sumar un 1.
//Obtengo la posicion donde terminan los hijos del nodo id
id = mapeo1select(id+1,tipo,0);
comprobar[id] = 1;
q.push(id);
for(int i = id-1; i >=0; i--){
    if(mapeodirecto1(i,tipo) == 0)break;
    if(mapeodirecto1(i,tipo)== 1){
        if(agregados[mapeo1rank(i+1,tipo,1)] == 0)sons.push_back(mapeo1rank(i+1,tipo,1));
        agregados[mapeo1rank(i+1,tipo,1)]++;
    }
    else {
        if(agregados[mapeodirecto2(mapeo1rank(i+1,tipo,2)-1,tipo)] == 0) sons.push_back(mapeodirecto2(mapeo1rank(i+1,tipo,2)-1,tipo));
        agregados[mapeodirecto2(mapeo1rank(i+1,tipo,2)-1,tipo)]++;
    }

}


while(!q.empty())q.pop();
vector <int> papas;
map<int,int> comprobar2;
//Obtengo su posicion en el wt que representa la secuencia B
// con esta linea mapeo de posicion en B a identificador en B
// para mapear de identificador a posicion se realiza un select(identificador,1)
id = idoriginal;
comprobar2[id] = 1;
q.push(id);
   if(id != 0) while(!q.empty()){
    id = q.front();
    q.pop();
    if(comprobar2[mapeo1rank(mapeo1select(id,tipo,1),tipo,0) ] == 0 ){
       if(agregados[mapeo1rank(mapeo1select(id,tipo,1),tipo,0)] == 0 ) sons.push_back(mapeo1rank(mapeo1select(id,tipo,1),tipo,0) );
       agregados[mapeo1rank(mapeo1select(id,tipo,1),tipo,0)]++;
        comprobar2[mapeo1rank(mapeo1select(id,tipo,1),tipo,0)]++;
    }
    int j = 1, limite = mapeo2rank(tam2(tipo),tipo,id);
    while(j <= limite){
        if(comprobar2[mapeo1rank(mapeo1select(mapeo2select(j,tipo,id)+1,tipo,2),tipo,0)] == 0 ){
            if(agregados[mapeo1rank(mapeo1select(mapeo2select(j,tipo,id)+1,tipo,2),tipo,0)] == 0)sons.push_back(mapeo1rank(mapeo1select(mapeo2select(j,tipo,id)+1,tipo,2),tipo,0));
            agregados[mapeo1rank(mapeo1select(mapeo2select(j,tipo,id)+1,tipo,2),tipo,0)]++;
            comprobar2[mapeo1rank(mapeo1select(mapeo2select(j,tipo,id)+1,tipo,2),tipo,0)]++;
        }
        j++;
    }
break;
    }

return sons;

}



vector <int>  vecino2(int id, int tipo){
 vector <int> sons;
if(limitenodos <= id)return sons;
map<int,int> comprobar;
queue<int> q;   
comprobar[0]++;
int original = id;


id = mapeo1select(id,tipo,1)+1;
//cout << "pase id" << endl;
comprobar[id] = 1;
q.push(id);
while(true){
    if(mapeodirecto1(id,tipo) == 0)break;
    else if(mapeodirecto1(id,tipo)== 1){
        if(comprobar[mapeo1rank(id+1,tipo,1)] == 0)sons.push_back(mapeo1rank(id+1,tipo,1));
        comprobar[mapeo1rank(id+1,tipo,1)]++; 
        id = mapeogetclose(id,tipo)+1;       
    }
    else{
        if(comprobar[mapeodirecto2(mapeo1rank(id+1,tipo,2)-1,tipo)] == 0) sons.push_back(mapeodirecto2(mapeo1rank(id+1,tipo,2)-1,tipo));
        comprobar[mapeodirecto2(mapeo1rank(id+1,tipo,2)-1,tipo)]++;
        id = id + 2;        
    }

}


int posact = mapeogetenclose(mapeo1selectbopen(original-1,tipo),tipo),nuevoB;
if(posact != tam3(tipo)){
nuevoB = mapeo1rankbopen(posact,tipo); // rank es exclusivo, en vez de restar 1 para simular indices desde 0 simplemente no le sumo 1
if(comprobar[nuevoB ] == 0 ){
    sons.push_back(nuevoB );
    comprobar[nuevoB]++;
}
}

int j = 1, limite = mapeo2rank(tam2(tipo),tipo,original);
while(j <= limite){
    int mapeoHtobp = mapeoselectbp2(mapeo2select(j,tipo,original),tipo);
    int mapeoEnclose= mapeogetenclose(mapeoHtobp,tipo);
    if(mapeoEnclose != tam3(tipo)){
    nuevoB = mapeo1rankbopen(mapeoEnclose,tipo);
    if(comprobar[nuevoB] == 0 ){
        sons.push_back(nuevoB);
        comprobar[nuevoB]++;
    }
    }
    j++;
}



return sons;
}








vector <pair <int,int > > hijos4(int id, int sumar){
 vector <pair<int,int> > sons;
map<int,int> comprobar;
bit_vector auxcon(bps.size(),0);
queue<int> q;
auxcon[id] = 1;
q.push(bopen_select(id+1));
while(!q.empty()){
id = q.front();
q.pop();
 int posclose = bps.find_close(id);
sons.push_back(make_pair(bopen_rank(id+1),bopen_rank(posclose+1))); 
int limitel = b2H_rank(bps.rank(id)+1), limiter = b2H_rank(bps.rank(posclose)+1);


vector <long unsigned int> maprangos; //=  restricted_unique_range_values(wt2,limitel,limiter-1,minwt2,maxwt2);
if(limiter > limitel) maprangos =  restricted_unique_range_values(wt2mapedS,limitel,limiter - 1,minwt2,maxwt2);

for(int i = 0; i < maprangos.size();i++){

int evaluar = maprangos[i];

if(auxcon[evaluar] == 0 &&  mapbarato[evaluar] == 0){
    auxcon[evaluar] = 1;
    q.push(evaluar);
}
else if(auxcon[evaluar] == 0  &&   mapbarato[evaluar] == 1){
    auxcon[evaluar] = 1;
    q.push(evaluar);
    sons.push_back(make_pair(bopen_rank(id+1),-1));
}

}
}
return sons;
}





vector <int> padres2(int id){
/*
Funcion para obtener todos los padres dado un nodo dada un id en B
Si quiero unicamente los padres directo simplemente quito el while (que itere unicamente 1 vez), dado que
un nodo en un arbol tiene unicamente un padre sera un retorno directo
*/
vector <int> papas;
if(id == 0)return papas;
map<int,int> comprobar;
queue<int> q;
//Obtengo su posicion en el wt que representa la secuencia B
// con esta linea mapeo de posicion en B a identificador en B
// para mapear de identificador a posicion se realiza un select(identificador,1)
comprobar[id] = 1;
q.push(bopen_select(id+1));
    while(!q.empty()){
    id = q.front();
    if(id == 0)break;
    q.pop();
int posact = bps.enclose(id);
    // cout << "padre: " << id << " " << bopen_select(id+1) << " " << posact << endl;
    int nuevoB = bopen_rank(posact); // rank es exclusivo, en vez de restar 1 para simular indices desde 0 simplemente no le sumo 1
    if(comprobar[nuevoB ] == 0 ){
        papas.push_back(nuevoB );
        q.push(bopen_select(nuevoB+1));
        comprobar[nuevoB]++;
    }
    //Limite se calcula asi: busco si para el nodo actualmente analizado existe alguna aparicion en H, de ser asi obtengo todas su apariciones
    int j = 1, limite = wt2mapedS.rank(wt2mapedS.size(),id);
    while(j <= limite){
int mapeoHtobp = wt1.select(wt2mapedS.select(j,id)+1,2);
        nuevoB = bopen_rank(bps.enclose(mapeoHtobp));
        if(comprobar[nuevoB] == 0 ){
            papas.push_back(nuevoB);
            q.push(bopen_select(nuevoB+1));
            comprobar[nuevoB]++;
        }
        j++;
    }
    }
    return papas;

}





vector <int> hijos2(int id){
 vector <int> sons;
map<int,int> comprobar;
queue<int> q;   
comprobar[id]++;
q.push(id);

while(!q.empty()){
id = q.front();
id = bopen_select(id+1);
int i = bps.rank(id)+1;
int posclose = bps.find_close(id);
int limitef = bps.rank(posclose);
q.pop();

for( ; i <= limitef; i++){
int posact = bps.select(i);
if(bp2[posact] == 1){
int nodoactual = b2_rank(posact);
if(comprobar[wt2[nodoactual]] == 0){
    comprobar[wt2[nodoactual]]++;
    sons.push_back(wt2[nodoactual]);
    q.push(wt2[nodoactual]);
}
}
else{
int nodoactual = bopen_rank(posact);
if(comprobar[nodoactual] == 0){
    comprobar[nodoactual]++;
    sons.push_back(nodoactual);
}

}

}

}

return sons;
}





//ID es identificador de un nodo en B, no su posicion
vector <int> hijos(int id){
 vector <int> sons;
map<int,int> comprobar;
queue<int> q;   
comprobar[id]++;
q.push(id);
while(!q.empty()){
id = q.front();
id = wt1.select(id+1,0);
q.pop();
for(int i = id-1; i >=0; i--){
    if(wt1[i] == 0)break;
    int nuevoB = wt1.rank(i+1,1);
    if(wt1[i] == 1 && comprobar[nuevoB] == 0){
        q.push(nuevoB);
        sons.push_back(nuevoB);
        comprobar[nuevoB]++;
        continue;
    }
    if(wt1[i] == 2)nuevoB =  wt2[wt1.rank(i+1,2)-1];
    if(wt1[i] == 2 && comprobar[nuevoB] == 0){
        q.push(nuevoB);
        sons.push_back(nuevoB);
        comprobar[nuevoB]++;
    }
}

}

return sons;
}



bool regla1(int g1, int g2, int total){
    vector <int> pg1;
    if(g1 != 0)pg1 = padres2(g1);
    for(int i = 0; i < pg1.size(); i++)if(pg1[i] == g2 )return true;
    return false;
}



bool regla1_2(int g1, int g2, int total){
vector <pair <int,int> > pg2;
pg2 = hijos4(g2,1);
for(int i = 0; i < pg2.size(); i++){
if(pg2[i].second == -1)pg2[i].second = pg2[i].first;
if(pg2[i].first <= g1 && pg2[i].second >= g1 )return true;
}
return false;
} 



/*
Obtengo padres y los mapeo todos a disjoint, luego marco los padres de g2 y compruebo si alguno de los padres mapeados de g1 tiene como vecino
alguno de los padres marcados previamente de g2, en caso de ser asi retorno true, false en caso contrario
*/
bool regla2(int g1, int g2, int total){
    vector <int> pg1;
    if(g1 != 0)pg1 = padres2(g1);
    vector <int >pg2;
    if(g2 != 0)pg2 = padres2(g2);
    for(int i = 0; i < pg1.size(); i++)pg1[i] = dids1[pg1[i]];
    for(int i = 0; i < pg2.size(); i++){
        pg2[i] = dids1[pg2[i]];
    }
    vector <bool> visitado(total,false);
    for(int i = 0 ; i < pg2.size();i++)visitado[pg2[i]] = true;
    for(int i = 0; i < pg1.size(); i++) {
        vector <int> pg1aux = vecino2(pg1[i],1);
        for(int j = 0; j < pg1aux.size();j++){
            if(visitado[pg1aux[j]] && pg1aux[j] != 0)  return true;
        }

    }
    return false;
}



bool regla4(int g1, int g2, int tam){

//Obtengo hijos y marco, si alguno de estos es marcado doble, retorno true, caso contrario false

    vector <int> pg1 = hijos2(g1),pg2 = hijos2(g2);

        vector <bool> visitado(tam,false);
    for(int i = 0; i < pg1.size(); i++) {
            visitado[pg1[i]] = true;
            //else return true;
    }
    for(int i = 0; i < pg2.size(); i++) {
            if(visitado[pg2[i]]) return true;
    }
    return false;
}




bool regla4_2(int g1, int g2, int tam){
vector < pair< int,int> > pg1 = hijos4(g1,1),pg2 = hijos4(g2,2);
vector <bool> visitado(tam,false);
vector <int> mapeadosg1,mapeadosg2;
pg1.push_back(make_pair(g1,g1+1));
pg2.push_back(make_pair(g2,g2+1));
vector <bool> comprobar(tam,false);

for(int i = 0; i < pg1.size();i++){
    if(pg1[i].second == -1 && !visitado[pg1[i].first]) {
        mapeadosg1.push_back(pg1[i].first);
        comprobar[pg1[i].first] = true;
        visitado[pg1[i].first] = true;
    }
    else if(pg1[i].second != -1) for(int j = pg1[i].first; j < pg1[i].second; j++){
        int auxmapeado =j; 
        if(!visitado[auxmapeado]){
            visitado[auxmapeado] = true;
            comprobar[auxmapeado] = true;
            mapeadosg1.push_back(auxmapeado);
        }
    } 
}
     vector <bool> visitado2(tam,false);

for(int i = 0; i < pg2.size();i++){
    if(pg2[i].second == -1 && !visitado2[pg2[i].first]) {
        mapeadosg2.push_back(pg2[i].first);
        if(comprobar[pg2[i].first])return true;
        comprobar[pg2[i].first] = true;
        visitado2[pg2[i].first] = true;
    }
    else if(pg2[i].second != -1) for(int j = pg2[i].first; j < pg2[i].second; j++){
        int auxmapeado = j;
        if(!visitado2[auxmapeado]){
            visitado2[auxmapeado] = true;
            if(comprobar[auxmapeado])return true;
            comprobar[auxmapeado] = true;
            mapeadosg2.push_back(auxmapeado);
        }
    } 
}
return false;
}

bool regla4_3(int g1, int g2, int tam){

    vector < pair< int,int> > pg1 = hijos4(g1,1),pg2 = hijos4(g2,2);



for(int i = 0; i < pg1.size(); i++){
    if(pg1[i].second == -1)pg1[i].second = pg1[i].first;
    for(int j = 0; j < pg2.size();j++){
        if(pg2[j].second == -1)pg2[j].second = pg2[j].first;
        if(pg1[i].first>= pg2[j].first && pg1[i].second <= pg2[j].second ){
            return true;
        }
        if(pg2[j].first>= pg1[i].first && pg2[j].second <= pg1[i].second ){
            return true;
        }
    }
}

        vector <bool> visitado(tam,false);
return false;
}



//Identico a regla2, solo que analizo hijos en vez de padres (si corrigo g2 correguir g5 igual)
bool regla5(int g1, int g2, int total){

vector < pair< int,int> > pg1 = hijos4(g1,1),pg2 = hijos4(g2,2);

vector <bool> visitado(total,false);
vector <int> mapeadosg1,mapeadosg2;
for(int i = 0; i < pg1.size();i++){
    if(pg1[i].second == -1 && !visitado[pg1[i].first]) {
        mapeadosg1.push_back(pg1[i].first);
        visitado[pg1[i].first] = true;
    }
    else if(pg1[i].second != -1) for(int j = pg1[i].first; j < pg1[i].second; j++){
        int auxmapeado = j;//permbarato[j];
        if(!visitado[auxmapeado]){
            visitado[auxmapeado] = true;
            mapeadosg1.push_back(auxmapeado);
        }
    } 
}
     vector <bool> visitado2(total,false);
for(int i = 0; i < pg2.size();i++){
    if(pg2[i].second == -1 && !visitado2[pg2[i].first]) {
        mapeadosg2.push_back(pg2[i].first);
        visitado2[pg2[i].first] = true;
    }
    else if(pg2[i].second != -1) for(int j = pg2[i].first; j < pg2[i].second; j++){
        int auxmapeado = j;//permbarato[j];
        if(!visitado2[auxmapeado]){
            visitado2[auxmapeado] = true;
            mapeadosg2.push_back(auxmapeado);
        }
    } 
}


    for(int i = 0; i < mapeadosg1.size(); i++)mapeadosg1[i] = dids3[mapeadosg1[i]];
    for(int i = 0; i < mapeadosg2.size(); i++)mapeadosg2[i] = dids3[mapeadosg2[i]];


    vector <bool> visitado3(total,false);
    for(int i = 0 ; i < mapeadosg2.size();i++)visitado3[mapeadosg2[i]] = true;
    for(int i = 0; i < mapeadosg1.size(); i++) {
        vector <int> pg1aux = vecino2(mapeadosg1[i],3);
        for(int j = 0; j < pg1aux.size();j++){
            if(visitado3[pg1aux[j]] && pg1aux[j] != 0)  return true;
        }

    }

    return false;
}



bool regla7_8(int g1, int g2, int total){

    vector <bool> visitados(total,false);
    vector <int> mapeadosg2;
    vector < pair< int,int> > pg2 = hijos4(g2,2);
    vector <int> pg1 = padres2(g1);
     vector <bool> visitado2(total,false);
for(int i = 0; i < pg2.size();i++){
    if(pg2[i].second == -1 && !visitado2[pg2[i].first]) {
        mapeadosg2.push_back(pg2[i].first);
        visitado2[pg2[i].first] = true;
    }
    else if(pg2[i].second != -1) for(int j = pg2[i].first; j < pg2[i].second; j++){
        int auxmapeado = j;//permbarato[j];
        if(!visitado2[auxmapeado]){
            visitado2[auxmapeado] = true;
            mapeadosg2.push_back(auxmapeado);
        }
    } 
}    
    for(int i = 0; i < pg1.size(); i++)pg1[i] = dids2[pg1[i]];
    for(int i = 0; i < mapeadosg2.size(); i++)mapeadosg2[i] = dids2[mapeadosg2[i]];
    for(int i = 0; i <  mapeadosg2.size();i++)visitados[mapeadosg2[i]] = true;
    for(int i = 0; i < pg1.size();i++){
        vector <int> pg1aux = vecino2(pg1[i],2);
        for(int j = 0; j < pg1aux.size();j++){
            if(visitados[pg1aux[j]] )  {
               return true;
            }
        }
    }
    return false;
}



bool reglasub(int a, int b){

if( regla1(a,b,mitotal +1) )return true;

return false;
}


bool regladis(int a, int b){

if( regla2(a,b,mitotal +1) )return true;

return false;
}


bool reglanotdis(int a, int b){

if( regla1(a,b,mitotal +1) )return true;

if( regla4_3(a,b,mitotal +1) )return true;

if( regla5(a,b,mitotal +1) )return true;

return false;
}

bool reglanotsub(int a, int b){

if( regla2(a,b,mitotal +1) )return true;

if( regla7_8(a,b,mitotal +1) )return true;

return false;
}

  };








}// end namespace sdsl
#endif
