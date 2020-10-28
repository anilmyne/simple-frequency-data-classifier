# include "algorithms.hpp"
# include <iostream>
# include <string>

#ifndef MARGINS_HPP_INCLUDED
#define MARGINS_HPP_INCLUDED

typedef unsigned int UINT;

template <class T>
class Com{
public:
    ~Com(){;}
    Com(std::string identifier, VECTOR_INT relation){ // memory & path
        identifier_ = identifier;
        relation_ = relation;
    }
    short getmode() const{ return this->mode;}
    void add(std::vector<T> items, bool setmgn=0, long mgn=-1);
    void add(T item, bool setmgn=0, long mgn=-1);
    void add(std::vector< std::vector<T> > items, bool setmgn=0, long mgn=-1);
    void getprobabilities();
    void marginalize();
    std::vector<T> getfacts() const{ return this->facts_0;}
    std::vector< std::vector<T> > getfacts_n() const{return this->facts_1;};
private:
    bool isupdate = 0, isrestricted = 0;
    std::string identifier_;
    unsigned long prevsum, maxmgn = (255 * 255 * 255) + 254;
    short mode, prec = 3;
    std::vector<int> relation_;
    std::vector<long> lastupdateindices;
    std::vector<T> mgnnames_0, facts_0;
    std::vector< std::vector<T> > mgnnames_1, facts_1;
    std::vector<long> mgns;
    std::vector<float> probs, probabilities;

};

template <class T>
void Com<T>::getprobabilities(){
    double sm = (double) std::accumulate(mgns.begin(), mgns.end(), 0);
    for (int i = 0; i < mgns.size(); i++){
        probabilities.push_back(atl::nround(mgns[i] / sm, prec));
    }
}

template <class T>
void Com<T>::add(T item, bool setmgn, long mgn){
    if (std::accumulate(mgns.begin(), mgns.end(), 0) == maxmgn){
        isupdate = 0;
        isrestricted = 1;
        return;
    }
    long pdiff = 0;
    mode = 0;
    lastupdateindices.clear();
    long indx = -1;
    for (long i = 0; i < mgnnames_0.size(); i++){
        if (mgnnames_0.at(i) == item){
            indx = i;
            break;
        }
    }
    if (indx != -1){
        lastupdateindices.push_back(indx);
        if (setmgn){
            pdiff = mgn - mgns[indx];
            mgns[indx] == mgn;
        }else{
            pdiff = 1;
            mgns[indx] += 1;
        }
    }else{
        lastupdateindices.push_back(mgnnames_0.size());
        mgnnames_0.push_back(item);
        if (setmgn){ pdiff += mgn; mgns.push_back(mgn); }
        else {mgns.push_back(1); pdiff = 1;}
    }
    prevsum = atl::sum(mgns) - pdiff;
    isupdate = 1;
}

template <class T>
void Com<T>::marginalize(){
    // finds STABLE probability
    if (!isupdate) return;
    isupdate = 0;
    long summ = atl::sum(mgns), usum = 0, nonupdatesum = 0;
    long num_upd;
    float urate, diff, avg;
    bool isvalid = 0;
    if (prevsum == 0) return;
    urate = 1 / (float)(summ * prevsum);
    num_upd = lastupdateindices.size();
    if (num_upd == 0) return;
    avg = (float) (summ / num_upd);
    //std::cout << "current: 1 / (" << summ  << "*" << prevsum << ") as " << urate << "|" << avg;
    isvalid = 1;
    VECTOR_FLOAT diffs, facts;
    long indx = -1;
    for (int i = 0; i < mgns.size(); i++){
        // use next indx match != 'in'
        for (int x=0; x < lastupdateindices.size(); x++){
            if (i == x){ indx = x; break;}
        }
        if (indx != -1){
            diff = urate * (avg - mgns[i]) * num_upd;
        }else{
            diff = urate * mgns[i] * num_upd;
        }
        diffs.push_back(diff);
        // print("diff:", diff, round(diff, self.prec))
        if (atl::nround(diff, prec) != 0){
            facts.clear();
            isvalid = 0;
            break;
        }
    }
    // print("isvalid:", isvalid)
    if (isvalid){
        Com::getprobabilities();
        float maxx = atl::max(probabilities);
        // print("MAXX:", maxx, self.probabilities)
        // r(maxx - pr[i]) == 0 ;
        for (int i = 0; i < probabilities.size(); i++){
            if (atl::nround(probabilities[i], prec) == maxx){
                if (mode == 0){
                    facts_0.push_back(mgnnames_0[i]);
                }else{
                    facts_1.push_back(mgnnames_1[i]);
                }
            }
        }
    }
}

template <class T>
class iCom{
public:
    ~iCom(){;}
    iCom(std::string identifier, std::vector<UINT> relation){ // memory & path
        identifier_ = identifier;
        relation_ = relation;
    }
    std::string id() const { return this->identifier_; }
    short getmode() const{ return this->mode;}
    void add(std::vector<T> items, bool setmgn=0, UINT mgn=0, bool modular=0);
    void add(T item, bool setmgn=0, UINT mgn=0);
    void add(std::vector< std::vector<T> > items, bool setmgn=0, UINT mgn=0, bool modular=0);
    void marginalize();
    void setbroadcast(bool b=1);
    bool isbroadcast();
    void write();
    void load();
    std::vector<T> getfacts() const{ return this->facts_0;}
    std::vector< std::vector<T> > getfacts_n() const{return this->facts_1;};
    std::vector<T> getnames() const{ return this->mgnnames_0; }
    std::vector<double> getipms() const{ return this->ipms; }  // t
    std::vector<UINT> getmgns() const{ return this->mgns; }
    std::vector<UINT> getinterrupts() const{ return this->interrupts; }
    void setnames(T name_){ this->mgnnames_0.push_back(name_); }
    void setipms(double ipm_){ this->ipms.push_back(ipm_); }  // t
    void setmgns(UINT mgn_){ this->mgns.push_back(mgn_); }
    void setinterrupts(UINT interrupt_){ this->interrupts.push_back(interrupt_); }
    bool testif (T item) const;
private:
    bool isupdate = 0, isrestricted = 0;
    std::string identifier_;
    UINT maxmgn = 255 * 255 * 255;
    short mode, prec = 3;
    std::vector<UINT> relation_;
    std::vector<T> mgnnames_0, facts_0;
    std::vector< std::vector<T> > mgnnames_1, facts_1;
    std::vector<UINT> mgns;
    std::vector<UINT> interrupts;
    std::vector<double> ipms;
    std::vector<double> ipmscache;;
};

template <class T>
void iCom<T>::add(T item, bool setmgn, UINT mgn){
    if (atl::sum(this->interrupts) >= this->maxmgn || atl::sum(this->mgns) >= this->maxmgn){
        this->facts_0 = this->mgnnames_0;
        this->facts_1 = this->mgnnames_1;
        this->isupdate = 0;
        this->isrestricted = 1;
        return;
    }
    this->isupdate = true;
    mode = 0;
    int indx = atl::binarySearch(this->mgnnames_0, item);
    if (indx != -1){
        if (setmgn){ this->mgns[indx] = mgn; }
        else{ this->mgns[indx] += 1; }
    }else{
        indx = atl::upperindexofappend(this->mgnnames_0, item);
        this->mgnnames_0.insert(this->mgnnames_0.begin() + indx, item);
        this->interrupts.insert(this->interrupts.begin() + indx, 0);
        this->ipms.insert(this->ipms.begin() + indx, 0.0);
        if (setmgn){ this->mgns.insert(this->mgns.begin() + indx, mgn); }
        else {this->mgns.insert(this->mgns.begin() + indx, 1); }
    }
    for (int i=0; i < this->mgns.size(); i++){
        if (i != indx){
            if (this->mgns[i] != 0){
                this->interrupts[i] += 1;
                this->ipms[i] += (1.0 / (double)(this->mgns[i]));
                this->mgns[i] = 0;
            }
        }
    }
}

template <class T>
void iCom<T>::add(std::vector<T> items, bool setmgn, UINT mgn, bool modular){
    if (atl::sum(this->interrupts) >= this->maxmgn || atl::sum(this->mgns) >= this->maxmgn){
        this->facts_0 = this->mgnnames_0;
        this->facts_1 = this->mgnnames_1;
        this->isupdate = 0;
        this->isrestricted = 1;
        return;
    }
    this->isupdate = true;
    mode = 0;
    std::vector<T> s_items;
    std::vector<int> lastupdateindices;
    int indx, lindx;
    for (int j=0; j < items.size(); j++){
        //std::cout << "\nSTART = "; atl::printvector<T>(this->mgnnames_0);
        indx = atl::binarySearch(this->mgnnames_0, items[j]);
        //std::cout << "\nindex of " << items[j] << " = " << indx;
        if (indx != -1){
            if (modular){
                lindx = atl::binarySearch(s_items, items[j]);
                if (lindx != -1){ continue; }
            }
            lindx = atl::upperindexofappend(s_items, items[j]);
            s_items.insert(s_items.begin()+lindx, items[j]);
            if (setmgn){ this->mgns[indx] = mgn; }
            else{ this->mgns[indx] += 1; }
        }else{
            if (modular){
                int lindx = atl::binarySearch(s_items, items[j]);
                if (lindx != -1){ continue; }
            }
            lindx = atl::upperindexofappend(s_items, items[j]);
            s_items.insert(s_items.begin()+lindx, items[j]);
            indx = atl::upperindexofappend(this->mgnnames_0, items[j]);
            this->mgnnames_0.insert(this->mgnnames_0.begin() + indx, items[j]);
            this->interrupts.insert(this->interrupts.begin() + indx, 0);
            this->ipms.insert(this->ipms.begin() + indx, 0.0);
            if (setmgn){ this->mgns.insert(this->mgns.begin() + indx, mgn); }
            else {this->mgns.insert(this->mgns.begin() + indx, 1); }
        }
    }
    for (int j=0; j < s_items.size(); j++){
        lastupdateindices.push_back(atl::binarySearch(this->mgnnames_0, s_items[j]));
    }
    for (int i=0; i < this->mgns.size(); i++){
        if (atl::binarySearch(lastupdateindices, i) != -1){ continue; }
        if (this->mgns[i] != 0){
            this->interrupts[i] += 1;
            this->ipms[i] += (1.0 / (double)(this->mgns[i]));
            this->mgns[i] = 0;
        }
    }
}

template <class T>
void iCom<T>::marginalize(){
    // finds nativeness : test with raw(no intr_div)
    std::vector<double> ipmcache;
    if (this->isupdate == false)
        return;
    for (int i=0; i < this->mgns.size(); i++){  // temporarily interrupt 0'd values
        double temp_ipm = this->ipms[i];
        if (this->mgns[i] != 0){
            temp_ipm += (1.0 / (double) mgns[i]);
            ipmcache.push_back(temp_ipm / (double)(this->interrupts[i] + 1)); // ADD HERE
        }else{
            ipmcache.push_back(temp_ipm / (double) this->interrupts[i]);  // avg
        }
    }
    double ipmcsum = 0, mid;
    for (int i=0; i<ipmcache.size(); i++){ ipmcsum += ipmcache[i]; }
    mid = ipmcsum / (double) ipmcache.size();
    //std::cout << "\nipmcache="; atl::printvector(ipmcache);
    //std::cout << "\nmid = " << mid << ":[" << ipmcsum << "]\n";
    for (int i=0;i < this->mgns.size(); i++){
        if (ipmcache[i] <= mid){  // less only?
            if (mode == 0){facts_0.push_back(this->mgnnames_0[i]);}
            else if(mode == 1){facts_1.push_back(this->mgnnames_1[i]);}
        }
    }
    this->ipmscache = ipmcache;
    isupdate = 0;
    //atl::printvector(mgns);
    //atl::printvector(interrupteds);
    //atl::printvector(interrupts);
    //atl::printvector(ipms);
    //std::cout << "<ipms\n";
}

template <class T>
bool iCom<T>::testif(T item) const{
    if (atl::binarySearch(this->facts_0, item) == -1){ return false;}
    return true;
}

template <class T>
void iCom<T>::load(){
    ;
}

template <class T>
void iCom<T>::write(){
    ;
}

template <class T>
class iCom_2{
public:
    ~iCom_2(){;}
    iCom_2(std::string identifier, VECTOR_INT relation){ // memory & path
        identifier_ = identifier;
        relation_ = relation;
    }
    short getmode() const{ return this->mode;}
    void add(std::vector<T> items, bool setmgn=0, int mgn=-1);
    void add(T item, bool setmgn=0, int mgn=-1);
    void add(std::vector< std::vector<T> > items, bool setmgn=0, int mgn=-1);
    void marginalize();
    std::vector<T> getfacts() const{ return this->facts_0;}
    std::vector< std::vector<T> > getfacts_n() const{return this->facts_1;};
private:
    bool isupdate = 0, isrestricted = 0;
    std::string identifier_;
    unsigned int maxmgn = (255 * 255 * 255) + 254;
    short mode, prec = 3;
    std::vector<int> relation_;
    std::vector<int> lastupdateindices;
    std::vector<T> mgnnames_0, facts_0;
    std::vector< std::vector<T> > mgnnames_1, facts_1;
    std::vector<int> mgns;
    std::vector<int> interrupts;
    std::vector<bool> interruptables;
};

template <class T>
void iCom_2<T>::add(T item, bool setmgn, int mgn){
    if (atl::sum(mgns) == maxmgn){
        facts_0 = mgnnames_0;
        facts_1 = mgnnames_1;
        isupdate = 0;
        isrestricted = 1;
        return;
    }
    isupdate = 1;
    mode = 0;
    lastupdateindices.clear();
    size_t indx = -1;
    for (size_t i=0; i < mgnnames_0.size(); i++){
        if (mgnnames_0[i] == item){
            indx = i;
            break;
        }
    }
    if (indx != -1){
        this->lastupdateindices.push_back(indx);
        if (setmgn){ this->mgns[indx] == mgn; }
        else{ mgns[indx] += 1; }
        this->interruptables[indx] = 1;
    }else{
        this->lastupdateindices.push_back(mgnnames_0.size());
        this->mgnnames_0.push_back(item);
        this->interrupts.push_back(0);
        this->interruptables.push_back(1);
        if (setmgn){ mgns.push_back(mgn); }
        else {mgns.push_back(1); }
    }
    for (int i=0; i < this->mgns.size(); i++){
            indx = -1;
        for (int x=0; x < this->lastupdateindices.size(); x++){
            if (lastupdateindices[x] == i){ indx = x; break; }
        }
        if (indx == -1 && this->interruptables[i] == 1){
            this->interrupts[i] += 1;
            this->interruptables[i] = 0;
        }
    }
}

template <class T>
void iCom_2<T>::marginalize(){
    // finds stability / randomity pattern
    std::vector<double> ipms;
    double dsum = 0;
    if (not this->isupdate) return;
    for (int i=0; i < this->mgns.size(); i++){  // temporarily interrupt 0'd values
        double divisor = 1;
        if (this->interrupts[i] != 0)
            divisor = (double) this->interrupts[i];
        dsum += divisor;
        ipms.push_back(divisor / (double) this->mgns[i]);
    }
    double mid = dsum / (double) (atl::sum(this->mgns));
    for (int i=0; i < this->mgns.size(); i++){
        if (ipms[i] < mid){
            if (mode == 0){ this->facts_0.push_back(this->mgnnames_0[i]); }
            else if(mode == 1){ this->facts_1.push_back(this->mgnnames_1[i]); }
        }
    }
    //atl::printvector(ipms);
    this->isupdate = 0;
}

#endif // MARGINS_HPP_INCLUDED
