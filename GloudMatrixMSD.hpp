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
wt_int< >      wt1,wtd1,wtnd1,wtns1,wt1p,wt2p,wt3p;
wt_int<> wt2,wtd2,wtnd2,wtns2,wt2maped,wt2mapedS;
    int_vector<> ids, ids2, ids3, ids4;
vector<int> original,original2,original3,original4;
vector <bool> inicial,inicial2,inicial3,inicial4;
vector <int> B,B2,B3,B4,H,H2,H3,H4,HM;
inv_perm_support<> perm, perm2, perm3, perm4;
int mitotal = 0;
bp_support_sada<256,32,rank_support_v<>,select_support_mcl<>> bps;
//bp2 bitvector marca 1 los 2, bopen marca 1 los 1
sd_vector<> bp2maped;
sd_vector<>::select_1_type bp2maped_select;
bit_vector Baux2, bp, bopen,bp2,hoja, undis,unndis,unnsub,reldis,relndis,relnsub,bp2H,mapbarato;
bit_vector::rank_1_type b_rank, b2_rank,bopen_rank, reldis_rank,relndis_rank,relnsub_rank,undis_rank,unndis_rank,unnsub_rank,b2H_rank ;
bit_vector::select_1_type b_select, b2_select,bopen_select,reldis_select,relndis_select,relnsub_select,undis_select,unndis_select,unnsub_select;
vector <int> wtbarato, permbarato;
int limitenodos,maxwt2,minwt2;

        void copy(const gbp& p) {
            original          = p.original;
            original2         = p.original2;
            original3         = p.original3;
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
            limitenodos     = p.limitenodos;
            maxwt2          = p.maxwt2;
            minwt2          = p.minwt2;

            HM  = p.HM;         
            Baux2 = p.Baux2;
            bp = p.bp; 
            bopen = p.bopen;
            bp2 = p.bp2;
            hoja = p.hoja; 

            undis = p.undis; 
            unndis = p.unndis; 
            unnsub = p.unnsub;
            reldis = p.reldis;
            relndis = p.relndis; 
            relnsub = p.relnsub;
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
    bit_vector auxA(B.size(),0),auxB(B.size(),0), auxC(B.size(),0),auxD(B.size(),0),auxE(B.size(),0),auxF(B.size(),0),auxI(B.size()/2,0),auxJ(B.size()+1,0),auxK(B.size(),0);
        int contadorabiertos = 0;
for(int i = 0 ; i < B.size(); i++){
    if(B[i] > 0)auxC[i] = 1;
    if(B[i] == 1)auxA[i] = 1;
    else if(B[i] == 2)auxB[i] = 1;
    if(B[i] == 2)auxD[i] = 1;
    if(B[i] == 1)auxE[i] = 1;
    if(B[i] == 2)auxI[contadorabiertos] = 1;
    if(B[i] == 1)contadorabiertos++;
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
sd_vector<> auxJJ(auxJ);
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

    vector <pair <int,int> > pares1;
    map <pair<int,int>,int> comprepe2;
    for(int i =  0; i < relaciones; i++){
        fscanf(fp2,"%d %d %d %d",&a,&b,&c,&d);
        if(comprepe2[make_pair(a+acumulado[b],c+acumulado[d])]== 1) continue;
        comprepe2[make_pair(a+acumulado[b],c+acumulado[d])]++;
        comprepe2[make_pair(c+acumulado[d],a+acumulado[b])]++;
        adjlist2[a+acumulado[b]].push_back(c+acumulado[d]);
        adjlist2[c+acumulado[d]].push_back(a+acumulado[b]);
        pares1.push_back( make_pair( mapOtoB(a+acumulado[b],0) , mapOtoB(c+acumulado[d],0)) );
        pares1.push_back( make_pair( mapOtoB(c+acumulado[d],0) , mapOtoB(a+acumulado[b],0)) );
    }


    sort(pares1.begin(), pares1.end());
    bit_vector auxexiste1(mitotal,0);
    vector<int> auxvun1;
    int actualrev = -1;
    for(int i = 0; i < pares1.size();i++){
        if(pares1[i].first != actualrev){
            actualrev = pares1[i].first;
            auxvun1.push_back(1);
            auxvun1.push_back(0);
            auxexiste1[actualrev] = 1;
        }
        else {
            auxvun1.push_back(0);
        }
    }
    auxvun1.push_back(1);
    reldis.swap(auxexiste1);
    bit_vector auxundis(auxvun1.size(),0);
    for(int i = 0 ; i < auxvun1.size();i++) if(auxvun1[i] == 1) auxundis[i] = 1;
    undis.swap(auxundis);



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


//    for(int i = 0; i < adjlist2[mitotal].size(); i++)cout << adjlist2[mitotal][i] << " ";
  //  cout << endl;
        queue <int> q;
    visitados2[mitotal] = true;
        vector<int> pdl(mitotal+1,-1);
    for(int i = 0; i < adjlist2[mitotal].size(); i++){
                int actual = adjlist2[mitotal][i];
               ids2[adjlist2[mitotal][i]] = posactual;
               posactual++;
        visitados2[actual] = true;
                q.push(actual);
                B2.push_back(1);
    }
    B2.push_back(0);
        while(!q.empty()){
            int frente = q.front();
   //         cout << "estoy en " << frente << endl;
            q.pop();
            for(int j = 0; j < adjlist2[frente].size();j++){
                if(adjlist2[frente][j] == pdl[frente])continue;
                if(!visitados2[adjlist2[frente][j]]){
                    visitados2[adjlist2[frente][j]] = true;
                    q.push(adjlist2[frente][j]);
                    pdl[adjlist2[frente][j]] = frente;
                    ids2[adjlist2[frente][j]] = posactual;
                    inicial2.push_back(false);
                    posactual++;
                    B2.push_back(1);
                    original.push_back(adjlist2[frente][j]);
                }
                else {
                    original3.push_back(adjlist2[frente][j]);
                    H2.push_back(ids2[adjlist2[frente][j]]);
            //        cout << "ID: " << ids[adjlist2[frente][j]] << endl;
                    B2.push_back(2);
                }

            }
            original2.push_back(frente);
            B2.push_back(0);

        }
    fclose(fp2);
    pdl.clear();


    FILE *fp3 = fopen(file3.c_str(), "r");
    if (!fp3) {
    fprintf(stderr, "Error opening file \"%s\".\n", file3.c_str());
    exit(EXIT_FAILURE);
    }
    fscanf(fp3,"%d",&relaciones);
    vector <pair <int,int> >pares2;
    vector <vector<int> > adjlist3(mitotal+1, vector<int>());
     map <pair<int,int>,int> comprepe;
    for(int i =  0; i < relaciones; i++){
        fscanf(fp3,"%d %d %d %d",&a,&b,&c,&d);
        if(comprepe[make_pair(a+acumulado[b],c+acumulado[d])]== 1) continue;
        comprepe[make_pair(a+acumulado[b],c+acumulado[d])]++;
        comprepe[make_pair(c+acumulado[d],a+acumulado[b])]++;
        adjlist3[a+acumulado[b]].push_back(c+acumulado[d]);
        adjlist3[c+acumulado[d]].push_back(a+acumulado[b]);
        pares2.push_back( make_pair( mapOtoB(a+acumulado[b],0) , mapOtoB(c+acumulado[d],0)) );
        pares2.push_back( make_pair( mapOtoB(c+acumulado[d],0) , mapOtoB(a+acumulado[b],0)) );
      //  cout << a+acumulado[b] << " "  <<c+acumulado[d] << endl;
    }



    sort(pares2.begin(), pares2.end());
    bit_vector auxexiste2(mitotal,0);
    vector<int> auxvun2;
    actualrev = -1;
    for(int i = 0; i < pares2.size();i++){
        if(pares2[i].first != actualrev){
            actualrev = pares2[i].first;
            auxvun2.push_back(1);
            auxvun2.push_back(0);
            auxexiste2[actualrev] = 1;
        }
        else {
            auxvun2.push_back(0);
        }
    }
    auxvun2.push_back(1);
    relnsub.swap(auxexiste2);
    bit_vector auxunsub(auxvun2.size(),0);
    for(int i = 0 ; i < auxvun2.size();i++) if(auxvun2[i] == 1) auxunsub[i] = 1;
    unnsub.swap(auxunsub);









    ids3.resize(mitotal+1);
    util::set_to_value(ids3,0);
    posactual = 1;
    visitados2.clear();
    unidos.clear();
    visitados2.resize(mitotal+1,false);
    unidos.resize(mitotal+1,false);

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


  //  for(int i = 0; i < adjlist3[mitotal].size(); i++)cout << adjlist3[mitotal][i] << " ";
    //cout << endl;
    visitados2[mitotal] = true;
    pdl.clear();
    pdl.resize(mitotal+1,-1);
    for(int i = 0; i < adjlist3[mitotal].size(); i++){
                int actual = adjlist3[mitotal][i];
               ids3[adjlist3[mitotal][i]] = posactual;
               posactual++;
        visitados2[actual] = true;
                q.push(actual);
                B3.push_back(1);
    }
    B3.push_back(0);
        while(!q.empty()){
            int frente = q.front();
      //      cout << "estoy en " << frente << endl;
            q.pop();
            for(int j = 0; j < adjlist3[frente].size();j++){
                if(adjlist3[frente][j] == pdl[frente])continue;
                if(!visitados2[adjlist3[frente][j]]){
                    visitados2[adjlist3[frente][j]] = true;
                    q.push(adjlist3[frente][j]);
                    pdl[adjlist3[frente][j]] = frente;
                    ids3[adjlist3[frente][j]] = posactual;
                    inicial3.push_back(false);
                    posactual++;
                    B3.push_back(1);
                    original.push_back(adjlist3[frente][j]);
                }
                else {
                    original3.push_back(adjlist3[frente][j]);
                    H3.push_back(ids3[adjlist3[frente][j]]);
            //        cout << "ID: " << ids[adjlist3[frente][j]] << endl;
                    B3.push_back(2);
                }

            }
            original3.push_back(frente);
            B3.push_back(0);

        }
    fclose(fp3);
    pdl.clear();

 FILE *fp4 = fopen(file4.c_str(), "r");
    if (!fp4) {
    fprintf(stderr, "Error opening file \"%s\".\n",file4.c_str());
    exit(EXIT_FAILURE);
    }
     fscanf(fp4,"%d",&relaciones);
    vector <vector<int> > adjlist4(mitotal+1, vector<int>());
    vector <pair <int,int> >pares3;
    map <pair<int,int>,int> mapita4;
    for(int i =  0; i < relaciones; i++){
        fscanf(fp4,"%d %d %d %d",&a,&b,&c,&d);
        if(mapita4[make_pair(a+acumulado[b],c+acumulado[d])]== 1) continue;
        mapita4[make_pair(a+acumulado[b],c+acumulado[d])]++;
        mapita4[make_pair(c+acumulado[d],a+acumulado[b])]++;
        adjlist4[a+acumulado[b]].push_back(c+acumulado[d]);
        adjlist4[c+acumulado[d]].push_back(a+acumulado[b]);
        pares3.push_back( make_pair( mapOtoB(a+acumulado[b],0) , mapOtoB(c+acumulado[d],0)) );
        pares3.push_back( make_pair( mapOtoB(c+acumulado[d],0) , mapOtoB(a+acumulado[b],0)) );

    }
    sort(pares3.begin(), pares3.end());
    bit_vector auxexiste3(mitotal,0);
    vector<int> auxvun3;
    actualrev = -1;
    for(int i = 0; i < pares3.size();i++){
        if(pares3[i].first != actualrev){
            actualrev = pares3[i].first;
            auxvun3.push_back(1);
            auxvun3.push_back(0);
            auxexiste3[actualrev] = 1;
        }
        else {
            auxvun3.push_back(0);
        }
    }
    auxvun3.push_back(1);
    relndis.swap(auxexiste3);
    bit_vector auxunndis(auxvun3.size(),0);
    for(int i = 0 ; i < auxvun3.size();i++) if(auxvun3[i] == 1) auxunndis[i] = 1;
    unndis.swap(auxunndis);


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
    pdl.clear();
    pdl.resize(mitotal+1,-1);
    for(int i = 0; i < adjlist4[mitotal].size(); i++){
                int actual = adjlist4[mitotal][i];
               ids4[adjlist4[mitotal][i]] = posactual;
               posactual++;
        visitados2[actual] = true;
                q.push(actual);
                B4.push_back(1);
    }
    B4.push_back(0);
        while(!q.empty()){
            int frente = q.front();
            q.pop();
            for(int j = 0; j < adjlist4[frente].size();j++){
                if(adjlist4[frente][j] == pdl[frente])continue;
                if(!visitados2[adjlist4[frente][j]]){
                    visitados2[adjlist4[frente][j]] = true;
                    q.push(adjlist4[frente][j]);
                    pdl[adjlist4[frente][j]] = frente;
                    ids4[adjlist4[frente][j]] = posactual;
                    inicial4.push_back(false);
                    posactual++;
                    B4.push_back(1);
                    original.push_back(adjlist4[frente][j]);
                }
                else {
                    original4.push_back(adjlist4[frente][j]);
                    H4.push_back(ids4[adjlist4[frente][j]]);
                    B4.push_back(2);
                }

            }
            original4.push_back(frente);
            B4.push_back(0);

        }
    fclose(fp4);
    pdl.clear();
util::init_support(reldis_rank, &reldis);
util::init_support(reldis_select, &reldis);
util::init_support(relndis_rank, &relndis);
util::init_support(relndis_select, &relndis);
util::init_support(relnsub_rank, &relnsub);
util::init_support(relnsub_select, &relnsub);


util::init_support(undis_rank, &undis);
util::init_support(undis_select, &undis);
util::init_support(unndis_rank, &unndis);
util::init_support(unndis_select, &unndis);
util::init_support(unnsub_rank, &unnsub);
util::init_support(unnsub_select, &unnsub);

util::init_support(reldis_rank, &reldis);
util::init_support(reldis_select, &reldis);
util::init_support(relndis_rank, &relndis);
util::init_support(relndis_select, &relndis);
util::init_support(relnsub_rank, &relnsub);
util::init_support(relnsub_select, &relnsub);



int_vector<> auxwt1(pares1.size()),auxwt2(pares2.size()),auxwt3(pares3.size());

for(int i = 0; i < pares1.size();i++) auxwt1[i] = pares1[i].second;

for(int i = 0; i < pares2.size();i++) auxwt2[i] = pares2[i].second;

for(int i = 0; i < pares3.size();i++) auxwt3[i] = pares3[i].second;

construct_im(wt1p,auxwt1);
construct_im(wt2p,auxwt2);
construct_im(wt3p,auxwt3);


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


//disjoint
construct_im(wtd1,vB2);
construct_im(wtd2,vH2);
//not subsumption
construct_im(wtns1,vB3);
construct_im(wtns2,vH3);
//not disjoint
construct_im(wtnd1,vB4);
construct_im(wtnd2,vH4);


//inv_perm_support<> auxperm1(&ids);
inv_perm_support<> auxperm2(&ids2);
inv_perm_support<> auxperm3(&ids3);
inv_perm_support<> auxperm4(&ids4);


//perm = auxperm1;
perm2 = auxperm2;
perm3 = auxperm3;
perm4 = auxperm4;

util::bit_compress(ids);
util::bit_compress(ids2);
util::bit_compress(ids3);
util::bit_compress(ids4);

for(int i = 0 ; i < ids.size(); i++) permbarato.push_back(perm[i]);
for(int i = 0; i < wt2.size();i++)wtbarato.push_back(wt2[i]);

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
            original3         = std::move(g.original3);
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

            undis = std::move(g.undis);
            unndis = std::move(g.unndis);
            unnsub = std::move(g.unnsub);
            reldis = std::move(g.reldis);
            relndis = std::move(g.relndis);
            relnsub = std::move(g.relnsub); 


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



      written_bytes += wt1p.serialize(out, child, "wt1p");
      written_bytes += wt1.serialize(out, child, "wt1");
      written_bytes += wt2p.serialize(out, child, "wt2p");
      written_bytes += wt3p.serialize(out, child, "wt3p");
      written_bytes += bps.serialize(out, child, "bps");


      written_bytes += reldis.serialize(out, child, "reldis");
      written_bytes += relndis.serialize(out, child, "relndis");
      written_bytes += relnsub.serialize(out, child, "relnsub");
      written_bytes += undis.serialize(out, child, "undis");
      written_bytes += unndis.serialize(out, child, "unndis");
      written_bytes += unnsub.serialize(out, child, "unnsub");
      written_bytes += undis_rank.serialize(out, child, "undis_rank");
      written_bytes += unndis_rank.serialize(out, child, "unndis_rank");
      written_bytes += unnsub_rank.serialize(out, child, "unnsub_rank");


      written_bytes += wt2mapedS.serialize(out, child, "wt2mapedS");

      written_bytes += undis_select.serialize(out, child, "undis_select");
      written_bytes += unndis_select.serialize(out, child, "unndis_select");
      written_bytes += unnsub_select.serialize(out, child, "unnsub_select");
      written_bytes += reldis_rank.serialize(out, child, "reldis_rank");
      written_bytes += relndis_rank.serialize(out, child, "relndis_rank");
      written_bytes += relnsub_rank.serialize(out, child, "relnsub_rank");
      written_bytes += reldis_select.serialize(out, child, "reldis_select");
      written_bytes += relndis_select.serialize(out, child, "relndis_select");
      written_bytes += relnsub_select.serialize(out, child, "relnsub_select");


      written_bytes += bopen.serialize(out, child, "bopen");
      written_bytes += bp2H.serialize(out, child, "bp2H");
      written_bytes += b2H_rank.serialize(out, child, "b2H_rank");
      written_bytes += bopen_rank.serialize(out, child, "bopen_rank");
      written_bytes += bopen_select.serialize(out, child, "bopen_select");


      written_bytes += perm.serialize(out, child, "perm");
      written_bytes += ids.serialize(out, child, "ids");


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
          
            Baux2.load(in);
            bp.load(in); 
            bopen.load(in);
            bp2.load(in);
            hoja.load(in); 

            undis.load(in); 
            unndis.load(in); 
            unnsub.load(in);
            reldis.load(in);
            relndis.load(in); 
            relnsub.load(in);

            mapbarato.load(in);
            bp2H.load(in);
            bp2maped.load(in);
        }


//RECORDAR QUE PARA TIPO = 2,3,4. SE PARTE DESDE 1 DADO QUE 0 ES NODO FICTICIO CREADO
pair <int,int> getrange(int id){
    int id2 = bopen_select(id+1);
    int posclose = bps.find_close(id2);
    return make_pair(id,bps.rank(posclose));
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
/*
Recibe par de ids de B que representan el inicio y el final y retorna
el mapeo al primer  y ultimo 1 que representa los ids mapeados en matriz
RETORNA ID DE B TRANSFORMADO A M(ASUMO QUE PARTEN EN 0 IGUAL QUE B)
*/
pair <int, int > mapBtoM(pair <int,int> rid, int tipo){
    int inicial, comprobarl, final, comprobarr;
    int frank,lrank,vtipo;
    if(tipo == 1)frank = reldis_rank(rid.first+1)-1;
    else if(tipo == 2) frank = relndis_rank(rid.first+1)-1;
    else frank = relnsub_rank(rid.first+1)-1;

    if(tipo == 1)vtipo = reldis[rid.first];
    else if(tipo == 2) vtipo = relndis[rid.first];
    else vtipo = relnsub[rid.first];

    inicial = frank;
    if(tipo == 1)final = reldis_rank(rid.second+1)-1; 
    else if(tipo == 2) final = relndis_rank(rid.second+1)-1;
    else final = relnsub_rank(rid.second+1)-1;
    if(inicial == final && vtipo == 0){
        inicial = -1;
        final = -1;
    }
    else if(inicial == -1){
        inicial = 0;
    }
 return make_pair(inicial,final);
}


/*
recibe id en M y retorna sus respectivas posiciones.
*/
pair <int, int> mapMIdtoMpos(pair <int,int> id, int tipo){
//asumo viene en indices desde el 0 al n-1
int inicial;
if(tipo == 1)inicial = undis_select(id.first+1);
else if(tipo == 2)inicial = unndis_select(id.first+1);
else inicial = unnsub_select(id.first+1);

if(inicial != 0) inicial -= (id.first);
int final;

if(tipo == 1)final =  undis_select(id.second+2) -  (id.second+2);
else if(tipo == 2)final = unndis_select(id.second+2) -  (id.second+2);
else final = unnsub_select(id.second+2) -  (id.second+2);
return make_pair(inicial,final);
}

 
//recibe un par correspondiente a posiciones en M y un segundo par correspondiente al rango de ID en B buscado
bool comprobarRango(pair <int,int>a,  pair<int,int> b , int tipo){
            if(tipo == 1){
        if(wt1p.range_search_2d(a.first,a.second,b.first,b.second,false).first > 0)return true;
    }
    if(tipo == 3)if(wt2p.range_search_2d(a.first,a.second,b.first,b.second,false).first > 0){
        return true;
    }        
    if(tipo == 2)if(wt3p.range_search_2d(a.first,a.second,b.first,b.second,false).first > 0)return true;
    return false;
}


//a y b son rangos a comprobar
bool comprobarRegla(vector <pair <int,int> >a,vector <pair <int,int> >b, int tipo ){
    for(int i = 0; i < a.size();i++){
        if(a[i].first == -1)continue;

        pair <int,int>auxfun = mapBtoM(a[i],tipo);
        if(auxfun.first == -1)continue;
       pair <int,int> pos = mapMIdtoMpos(auxfun,tipo );
        for(int j = 0; j < b.size();j++){
            if(b[j].first == -1 || b[j].second == -1)continue;

            pair<int,int>posb = b[j];
            if(posb.first < posb.second )b[j].second--;
            if( comprobarRango(pos,b[j],tipo)){
                return true;
            }
        }
    }
    return false;
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


//ID es identificador de un nodo en B, no su posicion
int mapeo1select(int id, int tipo, int x){
    if(tipo == 0) return wt1.select(id,x);
    if(tipo == 1) return wtd1.select(id,x);
    if(tipo == 2) return wtns1.select(id,x);
    if(tipo == 3) return wtnd1.select(id,x);
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
    //Limite se calcula asi: busco si para el nodo actualmente analizado existe alguna aparicion en H, de ser asi obtengo todas su apariciones
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




vector <pair <int,int > > hijos3(int id){
 vector <pair<int,int> > sons;
map<int,int> comprobar;
bit_vector auxcon(bps.size(),0);
queue<int> q;   
auxcon[id] = 1;
q.push(id);
while(!q.empty()){
id = q.front();
q.pop();
id = bopen_select(id+1);
int posclose = bps.find_close(id);
sons.push_back(make_pair(id,posclose));
int limitel = b2_rank(id+1), limiter = b2_rank(posclose+1); 

for(int i = limitel; i < limiter; i++){
int evaluar = wt2[i];
if(auxcon[evaluar] == 0 && hoja[i] == 0){
    auxcon[evaluar] = 1;
    q.push(evaluar);
}

}
//break;
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
//id = bopen_select(id+1);
int posclose = bps.find_close(id);
sons.push_back(make_pair(bopen_rank(id+1),bopen_rank(posclose+1)));
//int limitel = b2_rank(id+1), limiter = b2_rank(posclose+1); 
int limitel = b2H_rank(bps.rank(id)+1), limiter = b2H_rank(bps.rank(posclose)+1);

vector <long unsigned int> maprangos; //=  restricted_unique_range_values(wt2,limitel,limiter-1,minwt2,maxwt2);
if(limiter > limitel) maprangos =  restricted_unique_range_values(wt2maped,limitel,limiter - 1,minwt2,maxwt2);
//for(int i = limitel; i < limiter; i++){
for(int i = 0; i < maprangos.size();i++){
//int evaluar = wtbarato[i];

//int evaluar = maprangos[i];
int evaluar = bp2maped_select(maprangos[i]+1);
if(auxcon[evaluar] == 0 && /*hoja[i]*/ mapbarato[evaluar] == 0){
    auxcon[evaluar] = 1;
    q.push(evaluar);
}

else if(auxcon[evaluar] == 0  && /*hoja[i]*/  mapbarato[evaluar] == 1){
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
//     int posact =  bps.enclose(bopen_select(id+1));
    int posact = bps.enclose(id);
    // cout << "padre: " << id << " " << bopen_select(id+1) << " " << posact << endl;
    int nuevoB = bopen_rank(posact); // rank es exclusivo, en vez de restar 1 para simular indices desde 0 simplemente no le sumo 1
    if(comprobar[nuevoB ] == 0 ){
        papas.push_back(nuevoB );
   //     cout <<"padre1: " << nuevoB << endl;
        /*if(!inicial[nuevoB])*/q.push(bopen_select(nuevoB+1));
        comprobar[nuevoB]++;
    }
   // puts("sali");
    //Limite se calcula asi: busco si para el nodo actualmente analizado existe alguna aparicion en H, de ser asi obtengo todas su apariciones
//    int j = 1, limite = wt2.rank(wt2.size(),id);
        int j = 1, limite = wt2mapedS.rank(wt2mapedS.size(),id);
      //  puts("entre ");
    //    cout << j << " " << limite << " " <<id <<endl;
    while(j <= limite){
        //Antiguo
  //      int mapeoHtobp = b2_select(wt2.select(j,id)+1);
        //Alternativa 1
//        int mapeoHtobp = wt1.select(wt2.select(j,id)+1,2);
int mapeoHtobp = wt1.select(wt2mapedS.select(j,id)+1,2);
        //Alternativa 2
//        int mapeoHtobp = b_select(b2_select(2,wt2.select(j,id)+1)+1);//arreglar b_ actualmente b1 == b2
        nuevoB = bopen_rank(bps.enclose(mapeoHtobp));
        if(comprobar[nuevoB] == 0 ){
            papas.push_back(nuevoB);
        //            cout <<"padre2: " << nuevoB << endl;
            /*if(!inicial[nuevoB])*/q.push(bopen_select(nuevoB+1));
            comprobar[nuevoB]++;
        }
        j++;
    }
    //  puts("sali while ");
    }
   // puts("retornare ");
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
//cout << "entre con " << id << endl;
id = bopen_select(id+1);
int i = bps.rank(id)+1;
int posclose = bps.find_close(id);
int limitef = bps.rank(posclose);
//cout << i << " " << id << " " << posclose << " " << limitef << " " << bps.select(limitef) <<endl;
q.pop();
//int limitel = b2_rank(id), limiter = b2_rank(posclose); 

for( ; i <= limitef; i++){
int posact = bps.select(i);
if(bp2[posact] == 1){
int nodoactual = b2_rank(posact);
//cout << "nodo actual; " << nodoactual << " " << posact << endl;
if(comprobar[wt2[nodoactual]] == 0){
    comprobar[wt2[nodoactual]]++;
    sons.push_back(wt2[nodoactual]);
    q.push(wt2[nodoactual]);
}
}
else{
int nodoactual = bopen_rank(posact);
//cout << "entre if con " << nodoactual << endl;
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
//sumo 1 porque al ser un select e identificador partir de 0 debo sumar un 1.
//Obtengo la posicion donde terminan los hijos del nodo id
//id = wt1.select(id+1,0);
//cout << "consulto con id2 = " << id  << endl;
//cout << "padres de : " << id << endl;
while(!q.empty()){
id = q.front();
//cout << "entre con " << id << endl;
id = wt1.select(id+1,0);
//cout << "id: " << id << " " << wt1[id-1] << " " << wt1.size() <<endl;
q.pop();
for(int i = id-1; i >=0; i--){
  //  cout << wt1[i] << endl;
    if(wt1[i] == 0)break;
    int nuevoB = wt1.rank(i+1,1);
    if(wt1[i] == 1 && comprobar[nuevoB] == 0){
        q.push(nuevoB);
        sons.push_back(nuevoB);
        comprobar[nuevoB]++;
        continue;
        //cout << "marco " <<  wt1.rank(i+1,1) << endl;
    }
    if(wt1[i] == 2)nuevoB =  wt2[wt1.rank(i+1,2)-1];
    if(wt1[i] == 2 && comprobar[nuevoB] == 0){
       // cout << "entre al else " << endl;
        q.push(nuevoB);
        sons.push_back(nuevoB);
        comprobar[nuevoB]++;
         //       cout << "marco " <<  wt2[wt1.rank(i+1,2)-1] << endl;
    }
}
//cout << "termino for" << endl;

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



bool regla2_1(int g1, int g2, int tam){
//Obtengo hijos y marco, si alguno de estos es marcado doble, retorno true, caso contrario false
    vector <int> pg1;
    if(g1 != 0)pg1 = padres2(g1);
    vector <int> pg2;
    if(g2 != 0)pg2 = padres2(g2);
    vector <pair <int,int> > ppg1,ppg2;
    for(int i = 0; i < pg1.size();i++)ppg1.push_back(make_pair(pg1[i],pg1[i]));
    for(int i = 0; i < pg2.size();i++)ppg2.push_back(make_pair(pg2[i],pg2[i]));
    if(comprobarRegla(ppg1,ppg2,1))return true;
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
    for(int i = 0; i < pg1.size(); i++)pg1[i] = mapOtoB(mapBtoO(pg1[i],0),1);
    for(int i = 0; i < pg2.size(); i++){
        pg2[i] = mapOtoB(mapBtoO(pg2[i],0),1);
    }


    vector <bool> visitado(total,false);
    for(int i = 0 ; i < pg2.size();i++)visitado[pg2[i]] = true;
    for(int i = 0; i < pg1.size(); i++) {
        vector <int> pg1aux = vecinos(pg1[i],1);
        for(int j = 0; j < pg1aux.size();j++){
            if(visitado[pg1aux[j]] && pg1aux[j] != 0)  return true;
        }

    }
    return false;
}

bool regla4_3(int g1, int g2, int tam){
vector < pair< int,int> > pg1 = hijos4(g1,1),pg2 = hijos4(g2,2);
for(int i = 0; i < pg1.size(); i++){
    if(pg1[i].second == -1)pg1[i].second = pg1[i].first+1;
    for(int j = 0; j < pg2.size();j++){
        if(pg2[j].second == -1)pg2[j].second = pg2[j].first+1;
        if(pg1[i].first> pg2[j].first && pg1[i].second < pg2[j].second ){
      //      cout << g1 << " " << g2<<" "<<pg1[i].first<<" "<<pg1[i].second<<" "<<pg2[j].first <<" "<<pg2[j].second<<endl;
            return true;
        }
        if(pg2[j].first> pg1[i].first && pg2[j].second < pg1[i].second ){
        //    cout << g1 << " " << g2<<" "<<pg1[i].first<<" "<<pg1[i].second<<" "<<pg2[j].first <<" "<<pg2[j].second<<endl;
            return true;
        }
    }
}
return false;
}


bool regla4(int g1, int g2, int tam){

//Obtengo hijos y marco, si alguno de estos es marcado doble, retorno true, caso contrario false

    vector <int> pg1 = hijos2(g1),pg2 = hijos2(g2);
 
        vector <bool> visitado(tam,false);
        cout << pg1.size() << " " << pg2.size() << endl;
    for(int i = 0; i < pg1.size(); i++) {
            visitado[pg1[i]] = true;
            //else return true;
    }
    for(int i = 0; i < pg2.size(); i++) {
            if(visitado[pg2[i]]) return true;
    }
    return false;
}

bool regla4_1(int g1, int g2, int tam){
//cout << "entre con: " << g1 << " " << g2 << endl;
//Obtengo hijos y marco, si alguno de estos es marcado doble, retorno true, caso contrario false

    vector < pair< int,int> > pg1 = hijos4(g1,1),pg2 = hijos4(g2,2);
    for(int i = 0; i < pg1.size();i++)if(pg1[i].second == -1)pg1[i].second = pg1[i].first;
    for(int i = 0; i < pg2.size();i++)if(pg2[i].second == -1)pg2[i].second = pg2[i].first;
    //cout << pg1.size() <<" " << pg2.size() << endl;
    comprobarRegla(pg1,pg2,3);

   // cout << pg1.size() <<" " << pg2.size() << endl;
        vector <bool> visitado(tam,false);
return false;
/*
    for(int i = 0; i < pg1.size(); i++) {
            visitado[pg1[i]] = true;
            //else return true;
    }
    for(int i = 0; i < pg2.size(); i++) {
            if(visitado[pg2[i]]) return true;
    }
    return false;
    */
}


bool regla4_2(int g1, int g2, int tam){
vector < pair< int,int> > pg1 = hijos4(g1,1),pg2 = hijos4(g2,2);
vector <bool> visitado(tam,false);
vector <int> mapeadosg1,mapeadosg2;
vector <bool> comprobar(tam,false);

for(int i = 0; i < pg1.size();i++){
    if(pg1[i].second == -1 && !visitado[pg1[i].first]) {
        mapeadosg1.push_back(pg1[i].first);
        comprobar[pg1[i].first] = true;
        visitado[pg1[i].first] = true;
    }
    else if(pg1[i].second != -1) for(int j = pg1[i].first; j < pg1[i].second; j++){
        int auxmapeado =j; //permbarato[j];
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
        if(comprobar[pg2[i].first]){
           // cout << g1 << " " << g2<<endl;
            return true;
        }

        comprobar[pg2[i].first] = true;
        visitado2[pg2[i].first] = true;
    }
    else if(pg2[i].second != -1) for(int j = pg2[i].first; j < pg2[i].second; j++){
        int auxmapeado = j;//permbarato[j];
        if(!visitado2[auxmapeado]){
            visitado2[auxmapeado] = true;
            if(comprobar[auxmapeado]){
             //   cout << g1 << " " << g2<<endl;
                return true;
            }
            comprobar[auxmapeado] = true;
            mapeadosg2.push_back(auxmapeado);
        }
    } 
}
return false;
}


//Identico a regla2, solo que analizo hijos en vez de padres (si corrigo g2 correguir g5 igual)
bool regla5(int g1, int g2, int total){

    vector <int> pg1 = hijos2(g1),pg2 = hijos2(g2);
    for(int i = 0; i < pg1.size(); i++)pg1[i] = mapOtoB(mapBtoO(pg1[i],0),3);
    for(int i = 0; i < pg2.size(); i++)pg2[i] = mapOtoB(mapBtoO(pg2[i],0),3);

    vector <bool> visitado(total,false);
    for(int i = 0 ; i < pg2.size();i++)visitado[pg2[i]] = true;
    for(int i = 0; i < pg1.size(); i++) {
        vector <int> pg1aux = vecinos(pg1[i],3);
        for(int j = 0; j < pg1aux.size();j++){
            if(visitado[pg1aux[j]] && pg1aux[j] != 0)  return true;
//            else return true;
        }

    }
    return false;
}


bool regla5_1(int g1, int g2, int tam){
//Obtengo hijos y marco, si alguno de estos es marcado doble, retorno true, caso contrario false

    vector < pair< int,int> > pg1 = hijos4(g1,1),pg2 = hijos4(g2,2);
    for(int i = 0; i < pg1.size();i++)if(pg1[i].second == -1)pg1[i].second = pg1[i].first;
    for(int i = 0; i < pg2.size();i++)if(pg2[i].second == -1)pg2[i].second = pg2[i].first;
    if(comprobarRegla(pg1,pg2,2))return true;
return false;
}

bool regla7_8(int g1, int g2, int total){

    vector < pair< int,int> > pg2 = hijos4(g2,2);
    vector <int> pg1 = padres2(g1);

for(int i = 0; i < pg2.size();i++) if(pg2[i].first == pg2[i].second)pg2[i].first = -1;



    vector<pair <int,int> >ppadres;    
    for(int i = 0; i < pg2.size();i++)if(pg2[i].second == -1)pg2[i].second = pg2[i].first;
    for(int i = 0; i < pg1.size();i++) ppadres.push_back(make_pair(pg1[i],pg1[i]));
    if(comprobarRegla(ppadres,pg2,3))return true;
    return false;
}


/*
Obtengo todos los nodos con que g1 tiene la relacion directa de not subusmption, luego obtengo los padres de g2 y compruebo si alguno de estos
coincide con el set obtenido inicialmente.en caso de coincidir, se deriva la relacion.
*/
bool regla7(int g1, int g2, int total){
vector <int> pg2 = padres(g2);
return false;
    for(int i = 0; i < pg2.size(); i++)pg2[i] = mapOtoB(mapBtoO(pg2[i],0),2);
//                for(int i = 0; i < pg2.size(); i++) cout << pg2[i] << " ";
//    cout << endl;
    vector <bool> visitado(total,false);
    for(int i = 0 ; i < pg2.size();i++)visitado[pg2[i]] = true;
    vector <int> pg1aux = vecinos(mapOtoB(mapBtoO(g1,0),2),2);
    for(int j = 0; j < pg1aux.size();j++){
        if(visitado[pg1aux[j]] && pg1aux[j] != 0)  return true;
//            else return true;
        }
        return false;
}
/*
Obtengo los hijos de g1, y compruebo si alguno en el gloud tiene la relacion de not subsumption con g2
*/
bool regla8(int g1, int g2, int total){
vector <int> pg1 = hijos(g1);
//cout << g1 << " " << g2 << endl; 
  //      for(int i = 0; i < pg1.size(); i++) cout << pg1[i] << " ";
    //            cout << endl;
    for(int i = 0; i < pg1.size(); i++)pg1[i] = mapOtoB(mapBtoO(pg1[i],0),2);

   //             for(int i = 0; i < pg1.size(); i++) cout << pg1[i] << " ";
   // cout << endl;
        
    vector <bool> visitado(total,false);
    for(int i = 0 ; i < pg1.size();i++)visitado[pg1[i]] = true;
     //   cout << mapOtoB(mapBtoO(g2,0),2) << endl;
    vector <int> pg2aux = vecinos(mapOtoB(mapBtoO(g2,0),2),2);
    //            for(int i = 0; i < pg2aux.size(); i++) cout << pg2aux[i] << " ";
   // cout << endl;
    for(int j = 0; j < pg2aux.size();j++){
        if(visitado[pg2aux[j]] && pg2aux[j] != 0)  return true;
//            else return true;
        }
        return false;
}


/*
cambiar, la ida es que simule igual la regla de arriba, no sirve el range 2d query
*/

bool regla8_1(int g1, int g2, int tam){
//Obtengo hijos y marco, si alguno de estos es marcado doble, retorno true, caso contrario false

    vector < pair< int,int> > pg1 = hijos4(g1,1),pg2 = hijos4(g2,2);
    for(int i = 0; i < pg1.size();i++)if(pg1[i].second == -1)pg1[i].second = pg1[i].first;
    for(int i = 0; i < pg2.size();i++)if(pg2[i].second == -1)pg2[i].second = pg2[i].first;
    if(comprobarRegla(pg1,pg2,3))return true;

        vector <bool> visitado(tam,false);
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
