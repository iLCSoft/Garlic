#include <MCPhoton.hh>

ClassImp(MCPhoton)


MCPhoton::MCPhoton()
  : cosTheta(0),
    phi(0),
    Etot(0),
    E_GeV(0),
    RecEnRatio(0),
    nHits(0),
    zone(0), 
    interaction(0),
    rec(0),
    smallestDistToTrack(0),
    distToTrack(0),
    smallestDistToNextPhoton(0)
{
}

void MCPhoton::Streamer(TBuffer &R__b)
{
   // Stream an object of class MCPhoton.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> cosTheta;
      R__b >> phi;
      R__b >> Etot;
      R__b >> E_GeV;
      R__b >> RecEnRatio;
      R__b >> nHits;
      R__b >> zone;
      R__b >> interaction;
      R__b >> rec;
      R__b >> smallestDistToTrack;
      R__b >> distToTrack;
      R__b >> smallestDistToNextPhoton;
      R__b.CheckByteCount(R__s, R__c, MCPhoton::IsA());
   } else {
      R__c = R__b.WriteVersion(MCPhoton::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << cosTheta;
      R__b << phi;
      R__b << Etot;
      R__b << E_GeV;
      R__b << RecEnRatio;
      R__b << nHits;
      R__b << zone;
      R__b << interaction;
      R__b << rec;
      R__b << smallestDistToTrack;
      R__b << distToTrack;
      R__b << smallestDistToNextPhoton;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

