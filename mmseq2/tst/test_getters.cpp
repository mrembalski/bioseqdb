#include "gtest/gtest.h"
#include "mmseq2.h"

TEST(mockGetters, getSequence) {
    mock::init_mock_test({}, {"ACBDFE", "ADBECFA"}, 3);
    EXPECT_EQ(mock::get_indexes("NOT IMPORTANT", "ABC", 3), 0);

    mock::init_mock_test({}, {"ACB"}, 6);
    EXPECT_EQ(mock::get_indexes("BECAUSE THIS", "ABCDEF", 6), 0);

    mock::init_mock_test({}, {"ACBAA", "AAAA", "B", "DCAAFAB"}, 2);
    EXPECT_EQ(mock::get_indexes("GETTERS ARE", "AA", 2), 5);

    mock::init_mock_test({}, {"ABCDEFG", "AABCDEFA", "ABCDEFGABCDEFG", "DAGDAGDAG"}, 7);
    EXPECT_EQ(mock::get_indexes("MOCK FUNCTIONS", "ABCDEFG", 7), 3);
}

TEST(mockGetters, getIndexes) {
    uint64_t targedId;
    uint32_t position;

    mock::init_mock_test({}, {"ACBAA", "AAAA", "B", "DCAAFAB"}, 2);
    mock::get_ith_index(3, &targedId, &position, "AA", 2);
    EXPECT_EQ(targedId, 1);
    EXPECT_EQ(position, 2);

    mock::init_mock_test({}, {"ABCDEFG", "AABCDEFA", "ABCDEFGABCDEFG", "DAGDAGDAG"}, 7);
    mock::get_ith_index(2, &targedId, &position, "ABCDEFG", 7);
    EXPECT_EQ(targedId, 2);
    EXPECT_EQ(position, 7);
}

TEST(mockGetters, getIthIndex) {
    mock::init_mock_test({}, {"ABC", "D", "AB", "C"}, 1);
    EXPECT_EQ(*mock::get_sequence("TEST NAME", 2).get(), "AB");
}