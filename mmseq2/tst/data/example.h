#ifndef MMSEQ2_EXAMPLE_H
#define MMSEQ2_EXAMPLE_H

#include <string>
#include <algorithm>
#include <map>
#include <vector>

using VecStr = std::vector<std::string>;
using PairIdType = std::pair<std::string, std::string>;
using VecPairIdType = std::vector<PairIdType>;
using MapIdToSeqType = std::map<std::string, std::string>;

VecPairIdType getMatchesId() {
    return {
            {"A0A097J303", "D9IEU9"},
            {"A0A097J303", "C4N000"},
            {"A0A097J303", "A0A0K1LKA1"},
            {"A0A097J303", "Q7Y4P6"},
            {"A0A097J330", "A0A097J7L1"},
            {"A0A0K2QW23", "X1D7D6"},
            {"P30255", "P30257"},
            {"P59726", "C3VIT4"},
            {"P86485", "P62570"},
            {"Q1EN11", "F1AEL5"},
            {"Q1EN11", "Q17UY9"},
            {"Q1EN11", "L0P329"},
            {"Q1EN11", "P85882"},
            {"Q38L33", "S9BC69"},
            {"Q86QV4", "P0C892"},
            {"V5R462", "I6RTZ4"},
            {"W7V0Q8", "A0A059XRQ4"},
            {"W7V0Q8", "N9TVJ1"},
            {"W7V0Q8", "D3DIM0"},
            {"W7V0Q8", "O66982"},
            {"W7V0Q8", "A8UUM0"}
    };
}

MapIdToSeqType getQueryMapSequences() {
    return {
            {"A0A097J303", "MIKKILGYSLALATLLVALYYGVMFGLIQVVLFISDVIMALHSLVW"},
            {"A0A097J330", "MSSLWWCFVWLISIPLICLTFTFVMRLL"},
            {"A0A0K2QW23", "MIECEHLCMSMRGVRKPGAKTITSAVRGELRLPAARAEALSLIHGR"},
            {"P30255", "ENFAGGCTPGYQRTADGRCKATF"},
            {"P59726", "MWKKPAFIDLRLGLEVTLYISNR"},
            {"P86485", "GLVSSIGRALGGLLADVVKSKGQPA"},
            {"Q1EN11", "MAFLKKSLFLVLFLGLVSLSICDEEKRENEDEENQEDDEQSEMRRGLRSKIKEAAKTAGKMALGFVNDMAGEQ"},
            {"Q38L33", "MKNSKDILTNAIEEVSEKELMEVAGGKKGSGWFATITDDCPNSVFVCC"},
            {"Q86QV4", "DRDSCVDKSKCGKYGYYHQCDECCKKAGDRAGNCVYYKCKCNP"},
            {"V5R462", "MKYGVRYPKSGVHECMYGKRQAEQIWFLALRNGIAAVVVTDSGNGWEPS"}
    };
}

MapIdToSeqType getTargetMapSequences() {
    return {
            {"D9IEU9", "MIKKILGYSLALAALLVALYYGVMFGLIQVVLFISDVIMALHSLVW"},
            {"C4N000", "MIKKILGYSLALAALLVALYYGVIFGLIQVVLFISDVIMAIHSLVW"},
            {"A0A0K1LKA1", "MINKILGYSLALAALLVALYYGVMFGLIQVVLFISDVIMAIHSLVW"},
            {"Q7Y4P6", "MIKKILAGALGLLLLLTVLYYGVMFGLVQVVLFISDVIMVIHSLIW"},
            {"A0A097J7L1", "MSSLWWCFVWLISIPVICLTFTFVMRLL"},
            {"X1D7D6", "EHLCMSMRGVKKPNTLTVTSAVRGLFRENAASRAEVMALIKPSK"},
            {"P30257", "ENFVGGCTPGYQRTADGRCKPTF"},
            {"C3VIT4", "MWTKPSFEDLRLGLEVTLYISNR"},
            {"P62570", "GLVSSIGRALGGLLADVVKSKEQPA"},
            {"F1AEL5", "MFTLKKSLLLLFFIGVIKLSLCEEERNADEEKRRDDPDEMDVEVEKRLALERRDGWLRLFGLKPRRKH"},
            {"Q17UY9", "LVLFLGLVSLSICEEEKRETEEEENDQEEDDKSEEKRFLSLLPSIVSGAVSLAKKLG"},
            {"L0P329", "MDFLKKSLFLVVFLGLVSLSVCEEEKRESEEEKNEQEEDDREERSEEKRLLGMIPLAISAISALSKLG"},
            {"P85882", "MAFLKKSLFLVLFLGLVSLSICEEEKRETEEKENEQEDDDKSEEKRFLSLIPHAINAVSAIAKHFG"},
            {"S9BC69", "MKNSKDVLNNAIEEVSEKELMEVAGGKKGSGWFATITDDCPNSVFVCC"},
            {"P0C892", "DRDSCVDKSRCSKYGYYQECQDCCKKAGHNGGTCMFFKCKCA"},
            {"I6RTZ4", "MKWGVRYPISGVHECPFGKKQAEMIWFLAIRYGVRDPEVVVNHGDGVWRGAMK"},
            {"A0A059XRQ4", "MPKMKTKSALKKRIKITGTGKVLREQAYRSHLAQNKSTKQKRQARKSVQMHASDIKRFKGLF"},
            {"N9TVJ1", "MPKMKTKSALKKRIKITGTGKVLREQAYRSHLAQNKTTKQKRQARKSVQMHASDIKRFKGLF"},
            {"D3DIM0", "MAKVKMKSNRSAKKRFKITAKGKIKRWHAGGSHYNTKKAKDRKRRLRKPTLVNSGWEDKIRGLLKE"},
            {"O66982", "MAKVKMKTNRSAAKRFKVTAKGKIKRWKSGGAHYNTKKSSKRKRHLRKHTYVKDNMLKHVKALLKEF"},
            {"A8UUM0", "MAKVKMKTNRSAKKRFKVTATGKIKRWKGGGSHYNTKKDPKRKRRLRKATYVKENMEKHVRDLLRI"},
    };
}

VecStr getQueryIds() {
    auto mp = getQueryMapSequences();
    VecStr res;
    for (const auto &el : mp) {
        res.push_back(el.first);
    }
    return res;
}

VecStr getQuerySequences() {
    auto mp = getQueryMapSequences();
    VecStr res;
    for (const auto &el : mp) {
        res.push_back(el.second);
    }
    return res;
}

VecStr getTargetIds() {
    auto mp = getTargetMapSequences();
    VecStr res;
    for (const auto &el : mp) {
        res.push_back(el.first);
    }
    return res;
}

VecStr getTargetSequences() {
    auto mp = getTargetMapSequences();
    VecStr res;
    for (const auto &el : mp) {
        res.push_back(el.second);
    }
    return res;
}

#endif //MMSEQ2_EXAMPLE_H
