#include "catch.hpp"
#include "function.h"

using namespace std;
TEST_CASE("Basic Methods")
{
    HuberD* HD = new HuberD();
    CHECK(HD->get_v() == Approx(1.).epsilon(0.01));
    CHECK(HD->get_shift() == Approx(0.).epsilon(0.01));
    CHECK(HD->get_scale() == Approx(1.).epsilon(0.01));
}

TEST_CASE("Standard Distribution")
{
    HuberD* HD = new HuberD();
    CHECK(HD->Huber(0) == Approx(0.341).epsilon(0.01));
    CHECK(HD->Mksi_huber() == Approx(0).epsilon(0.01));
    CHECK(HD->Dksi_huber() == Approx(2.24).epsilon(0.01));
    CHECK(HD->asymmetry_huber() == Approx(0).epsilon(0.01));
    CHECK(HD->kurtosis_huber() == Approx(2.37).epsilon(0.01));
}

TEST_CASE("Primary; Shift and Scaled")
{
    HuberD* HD = new HuberD();
    HD->set_scale(2);
    HD->set_shift(2);
    CHECK(HD->Huber(0) == Approx(0.103).epsilon(0.01));
    CHECK(HD->Mksi_huber() == Approx(2).epsilon(0.01));
    CHECK(HD->Dksi_huber() == Approx(2.24).epsilon(0.01));
    CHECK(HD->asymmetry_huber() == Approx(0).epsilon(0.01));
    CHECK(HD->kurtosis_huber() == Approx(2.37).epsilon(0.01));
}

TEST_CASE("Mixture of Distributions")
{
    Mixture* MD = new Mixture();
    MD->get_component1()->set_scale(2);
    MD->get_component2()->set_scale(2);
    MD->get_component1()->set_shift(2);
    MD->get_component2()->set_shift(2);
    CHECK(MD->H_Mixture(0) == Approx(0.103).epsilon(0.01));
    CHECK(MD->Mksi_mixture() == Approx(2).epsilon(0.01));
    CHECK(MD->Dksi_Mixture() == Approx(2.24).epsilon(0.01));
    CHECK(MD->asymmetry_mixture() == Approx(0).epsilon(0.01));//?
    CHECK(MD->kurtosis_mixture() == Approx(1.77).epsilon(0.01));//?
}

TEST_CASE("Mixture Distribution - Expected value test")
{
    Mixture* MD = new Mixture();
    MD->get_component1()->set_scale(2);
    MD->get_component2()->set_scale(2);
    MD->get_component1()->set_shift(1);
    MD->get_component2()->set_shift(2);
    CHECK(MD->Mksi_mixture() == Approx(1.5).epsilon(0.01));
}

TEST_CASE("Mixture Distribution - Variance test")
{
    Mixture* MD = new Mixture();
    MD->get_component1()->set_scale(1);
    MD->get_component2()->set_scale(3);
    CHECK(MD->Dksi_Mixture() == Approx(2.24).epsilon(0.01));
}