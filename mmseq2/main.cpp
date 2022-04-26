#include "mmseq2.h"
#include "rpc/server.h"

int main()
{
    common::InputParams b;

    uint32_t port;
    std::cin >> port;

    rpc::server srv(port);

    srv.bind("mmseq2", [](common::InputParams inputParams)
             { return mmseq2::MMSeq2(inputParams); });

    srv.run();
    return 0;
}


//    runMMSeq2({"DDDDDDDDDCCGGGGGGGAA", "AAADDDDDDDCCGGGGGGGDD"},
//            {"DDDDDDDAAGGGGGGG"},
//            "blosum62", 7, 22, 0, 1, 4, 1, 2);

//    runMMSeq2({"AAADDDDDDDCCGGGGGGGDD"},
//            {"DDDDDDDAAGGGGGGG", "DDDDDDDDDCCGGGGGGGAA"},
//            "blosum62", 7, 22, 0, 1, 4, 1, 1);

//    runMMSeq2({"AACCTTGG", "ACTGACTGACTG", "TACTCAT"},
//            {"TACGGTAGCTTACTGA", "CTAGCTTACGATGCAAG", "CTTACAGCATACAGCATCGAT"},
//            "blosum62", 5, 10, 10, 1, 4, 1, 5);
