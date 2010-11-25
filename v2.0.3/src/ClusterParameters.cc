#include <ClusterParameters.hh>

ClassImp(ClusterParameters)


ClusterParameters::ClusterParameters()
  : ID(9999),
    theta(0),
    phi(0),
    Etot(0),
    Etot_g(0),
    Etot_g_noLC(0),
    Etot_g_phi(0),
    Etot_g_theta(0),
    Etot_g_pos(0),
    E_GeV(0),
    E_GeV_pre(0),
    E_GeV_pre_noPhi(0),
    E_GeV_noLC(0),
    E_GeV_noTheta(0),
    E_GeV_noPhi(0),
    E_GeV_en(0),
    E_GeV_hits(0),
    E_GeV_mix(0),
    E_GeV_opt(0),
    Es1(0),
    Es2(0),
    Es3(0),
    EnPs(0),
    En1odd(0),
    En1even(0),
    En2odd(0),
    En2even(0),
    NPs(0),
    N1odd(0),
    N1even(0),
    N2odd(0),
    N2even(0),
    start(0),
    end(0),
    depth(0),
    COGx(0),
    COGy(0),
    COGz(0),
    POSx(0),
    POSy(0),
    POSz(0),
    distToBiggest(0),
    smallestDistToBiggest(0),
    distFirstCell(0),
    distToTrack(0),
    smallestDistToTrack(0),
    E9C(0),
    E4C(0),
    E1C(0),
    nHits(0),
    nGhostHits(0),
    zone(0), 
    psHits(0),
    pl0Hits(0),
    hitDensity(0),
    enDensity(0),
    dirErr(0),
    seedDirErr(0),
    Eccentricity(0),
    Width(0),
    photonProb(0),
    photonProbMult(0),
    bckgrndProbMult(0),
    signalProbMult(0),
    surroundingLayers(0),
    particleID(0),
    isReal(0),
    MLP(0),
    finalCut(0),
    Volume(0),
    AssTo(0)
{
}

void ClusterParameters::Streamer(TBuffer &R__b)
{
   // Stream an object of class ClusterParameters.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> ID;
      R__b >> theta;
      R__b >> phi;
      R__b >> Etot;
      R__b >> Etot_g;
      R__b >> Etot_g_noLC;
      R__b >> Etot_g_phi;
      R__b >> Etot_g_theta;
      R__b >> Etot_g_pos;
      R__b >> E_GeV;
      R__b >> E_GeV_pre;
      R__b >> E_GeV_pre_noPhi;
      R__b >> E_GeV_noLC;
      R__b >> E_GeV_noTheta;
      R__b >> E_GeV_noPhi;
      R__b >> E_GeV_en;
      R__b >> E_GeV_hits;
      R__b >> E_GeV_mix;
      R__b >> E_GeV_opt;
      R__b >> Es1;
      R__b >> Es2;
      R__b >> Es3;
      R__b >> start;
      R__b >> end;
      R__b >> depth;
      R__b >> COGx;
      R__b >> COGy;
      R__b >> COGz;
      R__b >> POSx;
      R__b >> POSy;
      R__b >> POSz;
      R__b >> distToBiggest;
      R__b >> smallestDistToBiggest;
      R__b >> distFirstCell;
      R__b >> distToTrack;
      R__b >> smallestDistToTrack;
      R__b >> E9C;
      R__b >> E4C;
      R__b >> E1C;
      R__b >> EnPs;
      R__b >> En1odd;
      R__b >> En1even;
      R__b >> En2odd;
      R__b >> En2even;
      R__b >> NPs;
      R__b >> N1odd;
      R__b >> N1even;
      R__b >> N2odd;
      R__b >> N2even;
      R__b >> Chi2_long;
      R__b >> nHits;
      R__b >> nGhostHits;
      R__b >> zone;
      R__b >> psHits;
      R__b >> pl0Hits;
      R__b >> hitDensity;
      R__b >> enDensity;
      R__b >> dirErr;
      R__b >> seedDirErr;
      R__b >> Eccentricity;
      R__b >> Width;
      R__b >> photonProb;
      R__b >> photonProbMult;
      R__b >> bckgrndProbMult;
      R__b >> signalProbMult;
      R__b >> surroundingLayers;
      R__b >> particleID;
      R__b >> isReal;
      R__b >> MLP;
      R__b >> finalCut;
      R__b >> Volume;
      R__b >> AssTo;
      R__b.CheckByteCount(R__s, R__c, ClusterParameters::IsA());
   } else {
      R__c = R__b.WriteVersion(ClusterParameters::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << ID;
      R__b << theta;
      R__b << phi;
      R__b << Etot;
      R__b << Etot_g;
      R__b << Etot_g_noLC;
      R__b << Etot_g_phi;
      R__b << Etot_g_theta;
      R__b << Etot_g_pos;
      R__b << E_GeV;
      R__b << E_GeV_pre;
      R__b << E_GeV_pre_noPhi;
      R__b << E_GeV_noLC;
      R__b << E_GeV_noTheta;
      R__b << E_GeV_noPhi;
      R__b << E_GeV_en;
      R__b << E_GeV_hits;
      R__b << E_GeV_mix;
      R__b << E_GeV_opt;
      R__b << Es1;
      R__b << Es2;
      R__b << Es3;
      R__b << EnPs;
      R__b << En1odd;
      R__b << En1even;
      R__b << En2odd;
      R__b << En2even;
      R__b << NPs;
      R__b << N1odd;
      R__b << N1even;
      R__b << N2odd;
      R__b << N2even;
      R__b << start;
      R__b << end;
      R__b << depth;
      R__b << COGx;
      R__b << COGy;
      R__b << COGz;
      R__b << POSx;
      R__b << POSy;
      R__b << POSz;
      R__b << distToBiggest;
      R__b << smallestDistToBiggest;
      R__b << distFirstCell;
      R__b << distToTrack;
      R__b << smallestDistToTrack;
      R__b << E9C;
      R__b << E4C;
      R__b << E1C;
      R__b << Chi2_long;
      R__b << nHits;
      R__b << nGhostHits;
      R__b << zone;
      R__b << psHits;
      R__b << pl0Hits;
      R__b << hitDensity;
      R__b << enDensity;
      R__b << dirErr;
      R__b << seedDirErr;
      R__b << Eccentricity;
      R__b << Width;
      R__b << photonProb;
      R__b << photonProbMult;
      R__b << bckgrndProbMult;
      R__b << signalProbMult;
      R__b << surroundingLayers;
      R__b << particleID;
      R__b << isReal;
      R__b << MLP;
      R__b << finalCut;
      R__b << Volume;
      R__b << AssTo;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

