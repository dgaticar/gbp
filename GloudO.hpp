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




class gloud{
    public:
        typedef int_vector<>::size_type              size_type;
        typedef int_vector<>::value_type             value_type;
        typedef random_access_const_iterator<gloud> const_iterator;
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
wt_int<> wt2,wtd2,wtnd2,wtns2;
int_vector<> ids, ids2, ids3, ids4,wtbarato   , dids1,dids2,dids3;
vector<int> original,original2,original3,original4;
vector <bool> inicial,inicial2,inicial3,inicial4;
vector <int> B,B2,B3,B4,H,H2,H3,H4;
inv_perm_support<> perm, perm2, perm3, perm4;
int mitotal = 0;

        void copy(const gloud& p) {
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
        }


    public:

        //! Default constructor
        gloud() {};
//sep default = 1;
     	gloud(string file1, string file2, string file3, string file4) {
  int granularidades, aux,relaciones,a,b,c,d;
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
    for(int i = 0; i < totalnivel[granularidades-1]; i++){

        int actual = i + acumulado[granularidades-1];
        visitados[actual] = true;
        ids[actual] = posactual;
        inicial.push_back(true);
        posactual++;
        queue <int> q;
        q.push(actual);
        while(!q.empty()){
            int frente = q.front();
            q.pop();
            for(int j = 0; j < adjlist[frente].size();j++){
                if(!visitados[adjlist[frente][j]]){
                    visitados[adjlist[frente][j]] = true;
                    q.push(adjlist[frente][j]);
                    ids[adjlist[frente][j]] = posactual;
                    inicial.push_back(false);
                    posactual++;
                    B.push_back(1);
                    original.push_back(adjlist[frente][j]);
                }
                else {
                    original3.push_back(adjlist[frente][j]);
                    H.push_back(ids[adjlist[frente][j]]);
                    B.push_back(2);
                }

            }
            original2.push_back(frente);
            B.push_back(0);

        }
    }
    bit_vector auxA(B.size(),0),auxB(B.size(),0);
for(int i = 0 ; i < B.size(); i++){
    if(B[i] == 1)auxA[i] = 1;
    else if(B[i] == 2)auxB[i] = 1;
}


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
                    original2.push_back(adjlist2[frente][j]);
                }
                else {
                    original2.push_back(adjlist2[frente][j]);
                    H2.push_back(ids2[adjlist2[frente][j]]);
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
    vector <int> realceros;
    fscanf(fp3,"%d",&relaciones);
    vector <vector<int> > adjlist3(mitotal+1, vector<int>());
    for(int i =  0; i < relaciones; i++){
        fscanf(fp3,"%d %d %d %d",&a,&b,&c,&d);
            if(a+acumulado[b] == mitotal-1)realceros.push_back(c+acumulado[d]);
         if(c+acumulado[d]== mitotal-1)realceros.push_back(a+acumulado[b] );
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
    pdl.clear();
    pdl.resize(mitotal+1,-1); 
    queue <int> q2;
    for(int i = 0; i < adjlist3[mitotal].size(); i++){
                int actual = adjlist3[mitotal][i];
               ids3[adjlist3[mitotal][i]] = posactual;
               posactual++;
        visitados2[actual] = true;
                q2.push(actual);
                B3.push_back(1);
    }
    B3.push_back(0);
        while(!q2.empty()){
            int frente = q2.front();
            q2.pop();
            for(int j = 0; j < adjlist3[frente].size();j++){
                if(adjlist3[frente][j] == pdl[frente])continue;
                if(!visitados2[adjlist3[frente][j]]){
                    visitados2[adjlist3[frente][j]] = true;
                    q2.push(adjlist3[frente][j]);
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

    int_vector<> vB(B.size()),vB2(B2.size()),vB3(B3.size()) ,vB4(B4.size()); 
    int_vector <> vH(H.size()),vH2(H2.size()),vH3(H3.size()), vH4(H4.size()); 
    sd_vector <> prueba1(auxA);
    sd_vector <> prueba2(auxB);
    rrr_vector <63> prueba3(auxA);
    rrr_vector <63> prueba4(auxB);
    for(int i = 0 ; i < B.size(); i++)vB[i] = B[i];
    for(int i = 0 ; i < H.size(); i++)vH[i] = H[i];

    for(int i = 0 ; i < B2.size(); i++)vB2[i] = B2[i];
    for(int i = 0 ; i < H2.size(); i++)vH2[i] = H2[i];

    for(int i = 0 ; i < B3.size(); i++)vB3[i] = B3[i];
    for(int i = 0 ; i < H3.size(); i++)vH3[i] = H3[i];

    for(int i = 0 ; i < B4.size(); i++)vB4[i] = B4[i];
    for(int i = 0 ; i < H4.size(); i++)vH4[i] = H4[i];

construct_im(wt1,vB);
construct_im(wt2,vH);
//disjoint
construct_im(wtd1,vB2);
construct_im(wtd2,vH2);
//not subsumption
construct_im(wtns1,vB3);
construct_im(wtns2,vH3);
//not disjoint
construct_im(wtnd1,vB4);
construct_im(wtnd2,vH4);


inv_perm_support<> auxperm1(&ids);
inv_perm_support<> auxperm2(&ids2);
inv_perm_support<> auxperm3(&ids3);
inv_perm_support<> auxperm4(&ids4);


perm = auxperm1;
perm2 = auxperm2;
perm3 = auxperm3;
perm4 = auxperm4;

util::bit_compress(ids);
util::bit_compress(ids2);
util::bit_compress(ids3);
util::bit_compress(ids4);


wtbarato.resize(wt2.size());
for(int i = 0; i < wt2.size();i++)wtbarato[i] = wt2[i];
util::bit_compress(wtbarato);

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
        gloud(const gloud& g) {
            copy(g);
        }

        //! Copy constructor
        gloud(gloud&& g) {
            *this = std::move(g);
        }

        //! Assignment operator
        gloud& operator=(const gloud g) {
            if (this != &g) {
                copy(g);
            }
            return *this;
        }

        //! Assignment move operator
        gloud& operator=(gloud&& g) {
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
            }
            return *this;
        }

        //! Swap operator
        void swap(gloud& g) {/*
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



	  written_bytes += wtbarato.serialize(out, child, "wtbarato");
	  written_bytes += wt1.serialize(out, child, "wt1");
	  written_bytes += wtnd1.serialize(out, child, "wtnd1");
	  written_bytes += wtd1.serialize(out, child, "wtd1");

	  written_bytes += wtns2.serialize(out, child, "wtns2");
	  written_bytes += wtnd2.serialize(out, child, "wtnd2");
	  written_bytes += wt2.serialize(out, child, "wt2");
	  written_bytes += wtd2.serialize(out, child, "wtd2");


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
	    wtbarato.load(in);
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
        }


//RECORDAR QUE PARA TIPO = 2,3,4. SE PARTE DESDE 1 DADO QUE 0 ES NODO FICTICIO CREADO


/*
Dado un identificador en B, retornar el identificador mapeado a indices originales (indices de entrada )
tipo refiere a cual de los gloud, sinedo 1= contencion, 2 = disjoint, 3 = nocontencion, 4 = notdisjoint
*/
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
    else sons.push_back(wtbarato[wt1.rank(i+1,2)-1]);
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
    if(tipo == 0) return wtbarato[id];
    if(tipo == 1) return wtd2[id];
    if(tipo == 2) return wtns2[id];
    if(tipo == 3) return wtnd2[id];
    return -1;
}
 int tam2(int tipo){
    if(tipo == 0) return wtbarato.size();
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
        int agregar = mapeo1rank(i+1,tipo,1);
        if(agregados[agregar] == 0)sons.push_back(agregar);
        agregados[agregar]++;
    }
    else {
        int agregar = mapeodirecto2(mapeo1rank(i+1,tipo,2)-1,tipo);
        if(agregados[agregar] == 0) sons.push_back(agregar);
        agregados[agregar]++;
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
    id = id;
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


//ID es identificador de un nodo en B, no su posicion
vector <int> hijos(int id){
 vector <int> sons;
map<int,int> comprobar;
queue<int> q;   
comprobar[id]++;
q.push(id);
//sumo 1 porque al ser un select e identificador partir de 0 debo sumar un 1.
//Obtengo la posicion donde terminan los hijos del nodo id
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
    if(wt1[i] == 2)nuevoB =  wtbarato[wt1.rank(i+1,2)-1];
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
    if(g1 != 0)pg1 = padres(g1);
    for(int i = 0; i < pg1.size(); i++)if(pg1[i] == g2)return true;
    return false;
}

/*
Obtengo padres y los mapeo todos a disjoint, luego marco los padres de g2 y compruebo si alguno de los padres mapeados de g1 tiene como vecino
alguno de los padres marcados previamente de g2, en caso de ser asi retorno true, false en caso contrario
*/
bool regla2(int g1, int g2, int total){
    vector <int> pg1;
    if(g1 != 0)pg1 = padres(g1);
    vector <int >pg2;
    if(g2 != 0)pg2 = padres(g2);
    for(int i = 0; i < pg1.size(); i++)pg1[i] = dids1[pg1[i]];
    for(int i = 0; i < pg2.size(); i++){
        pg2[i] = dids1[pg2[i]];
    }
    vector <bool> visitado(total,false);
    for(int i = 0 ; i < pg2.size();i++)visitado[pg2[i]] = true;
    for(int i = 0; i < pg1.size(); i++) {
        vector <int> pg1aux = vecinos(pg1[i],1);
        for(int j = 0; j < pg1aux.size();j++){
            if(visitado[pg1aux[j]] && pg1aux[j] != 0)  return true; //SE DECUENTA 0 POR SER RAIZ FICTICIA
        }

    }
    return false;
}



bool regla4(int g1, int g2, int tam){

//Obtengo hijos y marco, si alguno de estos es marcado doble, retorno true, caso contrario false
    vector <int> pg1 = hijos(g1),pg2 = hijos(g2);
    pg1.push_back(g1);
    pg2.push_back(g2);
 //   cout << pg1.size() << " " << pg2.size() << endl;
        vector <bool> visitado(tam,false);
    for(int i = 0; i < pg1.size(); i++) {
            visitado[pg1[i]] = true;
    }
    for(int i = 0; i < pg2.size(); i++) {
            if(visitado[pg2[i]]) return true;
    }
    return false;
}


bool regla5(int g1, int g2, int total){

    vector <int> pg1 = hijos(g1),pg2 = hijos(g2);
    for(int i = 0; i < pg1.size(); i++)pg1[i] = dids3[pg1[i]];
    for(int i = 0; i < pg2.size(); i++)pg2[i] = dids3[pg2[i]];
    vector <bool> visitado(total,false);
    for(int i = 0 ; i < pg2.size();i++)visitado[pg2[i]] = true;
    for(int i = 0; i < pg1.size(); i++) {
        vector <int> pg1aux = vecinos(pg1[i],3);
        for(int j = 0; j < pg1aux.size();j++){
            if(visitado[pg1aux[j]] && pg1aux[j] != 0)  return true;
        }

    }
    return false;
}


bool regla7_8(int g1, int g2, int total){

    vector <bool> visitados(total,false);
    vector <int> pg1 = padres(g1), pg2 = hijos(g2);

    for(int i = 0; i < pg1.size(); i++)pg1[i] = dids2[pg1[i]];
    for(int i = 0; i < pg2.size(); i++)pg2[i] = dids2[pg2[i]];
    for(int i = 0; i <  pg2.size();i++)visitados[pg2[i]] = true;
    for(int i = 0; i < pg1.size();i++){

        vector <int> pg1aux = vecinos(pg1[i],2);

        for(int j = 0; j < pg1aux.size();j++){
            if(visitados[pg1aux[j]] )  return true;
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

if( regla4(a,b,mitotal +1) )return true;

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
