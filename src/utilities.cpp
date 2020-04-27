// parse CARD.json to individual homology files
// and check whether they still exist
//
// also load ARO to AMR family and AMR family to ARO datastructure
//

// ===========================================================================
// Classes
// ===========================================================================

class CARD {
    std::string CARDfp; 
    std::string CARD_version; 

  public:
    CARD (int,int);
    std::map<std::string, std::vector<std::string>> aroToAMRFamily() {return (width*height);}
    std::map<std::string, std::vector<std::string>> AMRFamilyToARO() {return (width*height);}
};

CARD::CARD(std::string CARDfp) {
  width = a;
  height = b;
}

// ===========================================================================
// Functions
// ===========================================================================


