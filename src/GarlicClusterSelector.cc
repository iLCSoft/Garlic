#include "GarlicClusterSelector.hh"

#include <iostream>
#include <fstream>
#include <assert.h>

using std::cout;
using std::endl;

void GarlicClusterSelector::readInCutValues(std::string cutFileName) {

  cout << "GarlicClusterSelector::readInCutValues: reading in cluster cut file " << cutFileName << endl;

  _varCuts.clear();
  std::string line;
  std::ifstream infile (cutFileName.c_str());

  int pdg;
  std::string varClass;
  std::string vname;
  std::string tag;
  std::string dummy;
  float percentile;
  std::vector < float > params;
  float par;

  std::pair < int, std::string > keyClass;
  std::pair < int, std::string > key;

  while ( std::getline( infile, line ) ) {
    if ( line[0]=='#' ) {
      continue;
    }
    if ( line.find( "END_OF_FILE" ) != std::string::npos ) break;

    std::stringstream sline(line);
    sline >> tag;
    if ( tag == "PDG" ) {
      sline >> pdg;
      keyClass.first = pdg;
      key.first = pdg;
    } else if ( tag == "CLASS" ) {
      sline >> varClass;
      keyClass.second = varClass;
    } else if ( tag == "VAR" ) {
      sline >> vname;

      if ( find( _varClasses[keyClass].begin(), _varClasses[keyClass].end(), vname )==_varClasses[keyClass].end() )
	_varClasses[keyClass].push_back( vname );

      key.second = vname;


    } else if ( tag == "EVTFRAC_PARS" ) {
      sline >> percentile;
      sline >> dummy;
      params.clear();
      while (sline) {
	par=987654321;
	sline >> par; 
	if ( par!=987654321 ) params.push_back(par);
      }

      if ( _varCuts.find(key)!=_varCuts.end() ) {
	_varCuts[key][percentile]=params;
      } else {
	std::map < float, std::vector < float > > temp;
	temp[percentile]=params;
	_varCuts[key]=temp;
      }

    } else {
      cout << "ERROR, unknown tag in cuts input file!" << tag << endl;
      assert(0);
    }

  }

  if (verbose) {
    cout << "var classes: " << endl;
    for ( std::map < std::pair < int, std::string > , std::vector < std::string > >::iterator itt=_varClasses.begin(); itt!=_varClasses.end(); itt++) {
      cout << itt->first.first << " " << itt->first.second << " : ";
      for (size_t i=0; i<itt->second.size(); i++) {
	cout << itt->second[i] << " ";
      }
      cout << endl;
    }
  }

  if (verbose) dump_cutvalues();

  return;
  }

std::map < std::string, float > GarlicClusterSelector::cluster_select(GarlicExtendedCluster* ecl) {

  // this one decides on the final classification, depending on the various cuts that have been passed or failed

  if (verbose) cout << "hello from cluster_select " << _pdg << endl;

  std::map < std::string , std::vector < int > > classSelections;

  for ( std::map < std::pair < int, std::string > , std::vector < std::string > >::iterator itt=_varClasses.begin(); itt!=_varClasses.end(); itt++) {
    if ( itt->first.first != _pdg ) continue; // wrong pdg
    classSelections[itt->first.second] = cluster_classSelection(ecl, itt->first.second);
  }

  std::map < std::string , float > classSelection;
  bool okTight(true);
  bool okLoose(true);
  bool okVLoose(true);

  int totNFail(0);

  for ( std::map < std::string , std::vector < int > >::iterator itt=classSelections.begin(); itt!=classSelections.end(); itt++) {
    int nTight(0);
    int nLoose(0);
    int nVLoose(0);
    int nFail(0);
    std::vector < int > sels = itt->second;
    int classScore(0);
    for (size_t i=0; i<sels.size(); i++) {
      classScore+=sels[i];
      switch (sels[i] ) {
      case TIGHT:
	nTight++;
	break;
      case LOOSE:
	nLoose++;
	break;
      case VERYLOOSE:
	nVLoose++;
	break;
      case REJECT:
	nFail++;
	break;
      default:
	cout << "GarlicClusterSelector::cluster_select ERROR: unknown selection code!" << sels[i] << endl;
      }
    }

    totNFail+=nFail;

    if ( float(nTight) / float(sels.size() ) < 0.66 )                okTight=false;
    if ( float(nTight+nLoose) / float(sels.size() ) < 0.66 )         okLoose=false;
    if ( float(nTight+nLoose+nVLoose) / float(sels.size() ) < 0.5 )  okVLoose=false;

    if ( verbose ) {
      cout << "class selection:" << itt->first << " T,L,VL,F: " << nTight << " " << nLoose << " " << nVLoose << " " << nFail << " / " << sels.size() << 
	" ; score:" << classScore << " okT, okL, okVL?" << okTight << " " << okLoose << " " << okVLoose << endl;
    }

    classSelection[itt->first]= float(classScore) / float(TIGHT*sels.size()); // fraction of max score

  }

  int sel(0);
  if      (okTight)  {
    if      ( totNFail==0 ) sel=TIGHT;
    else if ( totNFail <3 ) sel=LOOSE;
  } else if (okLoose)  {
    if      ( totNFail==0 ) sel=LOOSE;
    else if ( totNFail <3 ) sel=VERYLOOSE;
  } else if (okVLoose && totNFail<3) 
    sel=VERYLOOSE;

  if (verbose && sel>0)
    cout << "SELECTING cluster with flag " << sel << " ; total fails = " << totNFail << endl;


  classSelection["TOTAL"]=sel;

  return classSelection;
}


void GarlicClusterSelector::dump_cutvalues() {
  cout << "dumping cut values...:" << endl;
  for ( std::map < std::pair < int, std::string > , std::map < float, std::vector < float > > > ::iterator itt=_varCuts.begin(); itt!=_varCuts.end(); itt++) {
    cout << itt->first.first << " " << itt->first.second << endl;
    std::map < float, std::vector < float > > bb = itt->second;
    for ( std::map < float, std::vector < float > >::iterator jtt=bb.begin(); jtt!=bb.end(); jtt++) {
      cout << "---- " << jtt->first << " : ";
      for (size_t k=0; k<jtt->second.size(); k++) {
	cout << jtt->second[k] << ", ";
      }
      cout << endl;
    }
  }

  return;
}

std::vector < int > GarlicClusterSelector::cluster_classSelection(GarlicExtendedCluster* ecl, std::string className) {
  std::vector < int > selection;
  std::pair < int, std::string > gg(_pdg, className);
  if ( _varClasses.find( gg )==_varClasses.end() ) {
    cout << "GarlicClusterSelector::cluster_classSelection ERROR: cannot find entry for pdg " << _pdg << " class " << className << endl;
  } else {
    std::vector < std::string > class_varnames = _varClasses[gg];
    for (size_t i=0; i<class_varnames.size(); i++) {
      selection.push_back( cluster_select_variable ( ecl, class_varnames[i] ) );
    }
  } 
  return selection;
}


//-----------------------------------
// these are the detailed functions
//-----------------------------------

int GarlicClusterSelector::cluster_select_variable(GarlicExtendedCluster* ecl, std::string variable) {

  if ( _pdg!=22 && _pdg!=11 ) {
    cout << "GarlicClusterSelector::cluster_select_variable: WEIRD PDG " << _pdg << endl;
    assert (0);
  }

  std::pair < int, std::string > cutset(_pdg, variable);

  if ( _varCuts.find( cutset )==_varCuts.end() ) {
    cout << "GarlicClusterSelector::cluster_select_variable: ERROR could not find cuts for variable with name " << variable << " for PDG " << _pdg << endl;
    assert(0);
  }

  float value = ecl->getClusterProperty(variable);
    
  int sel(0);
  
  if ( _varCuts.find(cutset)==_varCuts.end() ) {
    cout << "GarlicClusterSelector::cluster_select_variable: ERROR, could not find cuts for pdg " << cutset.first << " , variable " << cutset.second << endl;
    assert (0);
  } 

  std::map < float, std::vector < float > > parameters = _varCuts[cutset];

  float perc[4]={0.01, 0.05, 0.95, 0.99}; // percentiles
  float cut[4]={0};

  

  float log_e(0);
  if ( ecl->getAssociatedTrack() ) {
    log_e = log10 ( ecl->getAssociatedTrack()->getTotalMomentum() );
  } else {
    log_e = log10 ( ecl->getEnergy() );
  }


  if (verbose) {
    cout << "hi from cluster_select_variable: pdg " << _pdg << " var " << variable << " energy, logE = " << ecl->getEnergy() << " " << log_e << endl;
    //    dump_cutvalues();
  }

  for (int ip=0; ip<4; ip++) {
    if ( parameters.find( perc[ip] ) == parameters.end() ) {
      cout << "GarlicClusterSelector::cluster_select_variable: ERROR cannot find parameters for percentile " << perc[ip] << endl;
      assert(0);
    }
    std::vector < float > cc = parameters[ perc[ip] ];
    for ( size_t i=0; i<cc.size(); i++) {
      cut[ip] += cc[i] * pow( log_e, i );
    }
  }

  if      ( value > cut[1] && value < cut[2] ) sel=TIGHT; // 2->98 %
  else if ( value > cut[0] && value < cut[3] ) sel=LOOSE; // 1->99 %
  else if ( value > ( cut[0] - 3*fabs(cut[0]-cut[1]) ) && 
	    value < ( cut[3] + 3*fabs(cut[3]-cut[2]) ) ) sel=VERYLOOSE;

  if (verbose) {
    cout << "value: " << value << " cuts: " << cut[0] << " " << cut[1] << " " << cut[2] << " " << cut[3] << " : sel = " << sel << endl;
  }

  return sel;
}


