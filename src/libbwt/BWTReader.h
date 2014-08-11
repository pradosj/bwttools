//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// BWTWriterBinary - Read a run length encoded binary BWT file from disk
//
#ifndef BWTREADERBINARY_H
#define BWTREADERBINARY_H

#include "Occurrence.h"
#include "EncodedString.h"
#include "BWTReader.h"
#include "bwt.h"


enum BWIOStage {IOS_NONE,IOS_HEADER,IOS_BWSTR,IOS_PC,IOS_OCC,IOS_DONE};
enum BWFlag {BWF_NOFMI = 0,BWF_HASFMI};

const uint16_t RLBWT_FILE_MAGIC = 0xCACA;
const uint16_t BWT_FILE_MAGIC = 0xEFEF;

class bwt;


/*!\class BWTReaderBinary
 * \brief Read a run length encoded BWT file from disk
 */
class BWTReaderBinary {
    public:
        BWTReaderBinary(const std::string& filename);
        ~BWTReaderBinary();

        //
        void read(bwt* pRLBWT);

        void readHeader(size_t& num_strings, size_t& num_symbols, BWFlag& flag);
        char readBWChar();
        void readRuns(std::vector<RLUnit>& out, size_t numRuns);

    private:
        std::istream* m_pReader;
        BWIOStage m_stage;
        RLUnit m_currRun;
        size_t m_numRunsOnDisk;
        size_t m_numRunsRead;
};

#endif
