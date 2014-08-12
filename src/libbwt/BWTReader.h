#ifndef BWTREADERBINARY_H
#define BWTREADERBINARY_H

#include <cstdint>
#include <string>
#include <vector>
#include "BWTReader.h"
#include "rlestr.h"
#include "bwt.h"


enum BWIOStage {IOS_NONE,IOS_HEADER,IOS_BWSTR,IOS_PC,IOS_OCC,IOS_DONE};
enum BWFlag {BWF_NOFMI = 0,BWF_HASFMI};

const uint16_t RLBWT_FILE_MAGIC = 0xCACA;
const uint16_t BWT_FILE_MAGIC = 0xEFEF;


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
