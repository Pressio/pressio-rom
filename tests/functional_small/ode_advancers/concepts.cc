
#include <gtest/gtest.h>
#include "pressio/ode_concepts.hpp"

struct MyState {};
using IV = double;

using namespace pressio::ode;

namespace{
  struct Obs1{
    void operator()(StepCount /*unused*/,
		    double /*unused*/,
		    const MyState & /*unused*/) const{}
  };

  struct Obs2{
    // wrong order
    void operator()(StepCount /*unused*/,
		    const MyState & /*unused*/,
		    double /*unused*/) const{}
  };

  struct Obs3{
    // wrong order
    void operator()(const MyState & /*unused*/,
		    StepCount /*unused*/,
		    double /*unused*/) const{}
  };
} //end anonym namespace

TEST(ode, concepts_state_observer)
{
  using namespace pressio::ode;

  static_assert( StateObserver<Obs1, double, MyState>::value, "");
  static_assert( !StateObserver<Obs2, double, MyState>::value, "");
  static_assert( !StateObserver<Obs3, double, MyState>::value, "");
}

// --------------------------------------------

#define CHECK_TRAIT1_TRUE(T_)  do { \
  static_assert(StepSizePolicy<T_, IV>::value, #T_ " should be valid (StepSizePolicy)"); \
  EXPECT_TRUE((StepSizePolicy<T_, IV>::value)) << #T_; \
} while(0)

#define CHECK_TRAIT1_FALSE(T_) do { \
  static_assert(!StepSizePolicy<T_, IV>::value, #T_ " should be invalid (StepSizePolicy)"); \
  EXPECT_FALSE((StepSizePolicy<T_, IV>::value)) << #T_; \
} while(0)

#define CHECK_TRAIT2_TRUE(T_)  do { \
  static_assert(StepSizePolicyWithReductionScheme<T_, IV>::value, #T_ " should be valid (WithReduction)"); \
  EXPECT_TRUE((StepSizePolicyWithReductionScheme<T_, IV>::value)) << #T_; \
} while(0)

#define CHECK_TRAIT2_FALSE(T_) do { \
  static_assert(!StepSizePolicyWithReductionScheme<T_, IV>::value, #T_ " should be invalid (WithReduction)"); \
  EXPECT_FALSE((StepSizePolicyWithReductionScheme<T_, IV>::value)) << #T_; \
} while(0)

// ====================== StepSizePolicy ======================
// Valid
struct SSP_ByRefVoidConst {
  void operator()(StepCount, StepStartAt<IV>, StepSize<IV>&) const {}
};
struct SSP_ByRefVoidConstNoexcept {
  void operator()(StepCount, StepStartAt<IV>, StepSize<IV>&) const noexcept {}
};
struct SSP_TemplatedExactlyRef {
  template<class U,
           std::enable_if_t<std::is_same<U, StepSize<IV>>::value, int> = 0>
  void operator()(StepCount, StepStartAt<IV>, U&) const {}
};

// Invalid
struct SSP_ByValue {
  void operator()(StepCount, StepStartAt<IV>, StepSize<IV>) const {}
};
struct SSP_ByConstRef {
  void operator()(StepCount, StepStartAt<IV>, const StepSize<IV>&) const {}
};
struct SSP_ByRvalueRef {
  void operator()(StepCount, StepStartAt<IV>, StepSize<IV>&&) const {}
};
struct SSP_ByRefNonVoid {
  int operator()(StepCount, StepStartAt<IV>, StepSize<IV>&) const { return 0; }
};
struct SSP_ByRefVoidNonConstOp {
  void operator()(StepCount, StepStartAt<IV>, StepSize<IV>&) {}
};
struct SSP_WrongIndVar {
  void operator()(StepCount, StepStartAt<float>, StepSize<IV>&) const {}
};
struct SSP_WrongStepSizeElem {
  void operator()(StepCount, StepStartAt<IV>, StepSize<float>&) const {}
};
struct SSP_ForwardingParam {
  template<class U>
  void operator()(StepCount, StepStartAt<IV>, U&&) const {}
};
struct SSP_OverloadSet {
  void operator()(StepCount, StepStartAt<IV>, StepSize<IV>&) const {}
  void operator()(StepCount, StepStartAt<IV>, const StepSize<IV>&) const {}
  void operator()(StepCount, StepStartAt<IV>, StepSize<IV>) const {}
};
struct SSP_NoCall {};

// Single TEST collecting all SSP cases
TEST(StepSizePolicyTrait, AllCases) {
  CHECK_TRAIT1_TRUE(SSP_ByRefVoidConst);
  CHECK_TRAIT1_TRUE(SSP_ByRefVoidConstNoexcept);
  CHECK_TRAIT1_TRUE(SSP_TemplatedExactlyRef);

  CHECK_TRAIT1_FALSE(SSP_ByValue);
  CHECK_TRAIT1_FALSE(SSP_ByConstRef);
  CHECK_TRAIT1_FALSE(SSP_ByRvalueRef);
  CHECK_TRAIT1_FALSE(SSP_ByRefNonVoid);
  CHECK_TRAIT1_FALSE(SSP_ByRefVoidNonConstOp);
  CHECK_TRAIT1_FALSE(SSP_WrongIndVar);
  CHECK_TRAIT1_FALSE(SSP_WrongStepSizeElem);
  CHECK_TRAIT1_FALSE(SSP_ForwardingParam);
  CHECK_TRAIT1_FALSE(SSP_OverloadSet);
  CHECK_TRAIT1_FALSE(SSP_NoCall);
}

// ====================== StepSizePolicyWithReductionScheme ======================
// Valid
struct SRP_ByRefVoidConst {
  void operator()(StepCount,
                  StepStartAt<IV>,
                  StepSize<IV>&,
                  StepSizeMinAllowedValue<IV>&,
                  StepSizeScalingFactor<IV>&) const {}
};
struct SRP_ByRefVoidConstNoexcept {
  void operator()(StepCount,
                  StepStartAt<IV>,
                  StepSize<IV>&,
                  StepSizeMinAllowedValue<IV>&,
                  StepSizeScalingFactor<IV>&) const noexcept {}
};
struct SRP_TemplatedExactlyRef {
  template<class A, class B, class C,
           std::enable_if_t<
             std::is_same<A, StepSize<IV>>::value &&
             std::is_same<B, StepSizeMinAllowedValue<IV>>::value &&
             std::is_same<C, StepSizeScalingFactor<IV>>::value,
             int> = 0>
  void operator()(StepCount, StepStartAt<IV>, A&, B&, C&) const {}
};

// Invalid: by-value (each param)
struct SRP_ByValue_dt {
  void operator()(StepCount, StepStartAt<IV>,
                  StepSize<IV>,
                  StepSizeMinAllowedValue<IV>&,
                  StepSizeScalingFactor<IV>&) const {}
};
struct SRP_ByValue_min {
  void operator()(StepCount, StepStartAt<IV>,
                  StepSize<IV>&,
                  StepSizeMinAllowedValue<IV>,
                  StepSizeScalingFactor<IV>&) const {}
};
struct SRP_ByValue_scale {
  void operator()(StepCount, StepStartAt<IV>,
                  StepSize<IV>&,
                  StepSizeMinAllowedValue<IV>&,
                  StepSizeScalingFactor<IV>) const {}
};

// Invalid: const& (each param)
struct SRP_ByConstRef_dt {
  void operator()(StepCount, StepStartAt<IV>,
                  const StepSize<IV>&,
                  StepSizeMinAllowedValue<IV>&,
                  StepSizeScalingFactor<IV>&) const {}
};
struct SRP_ByConstRef_min {
  void operator()(StepCount, StepStartAt<IV>,
                  StepSize<IV>&,
                  const StepSizeMinAllowedValue<IV>&,
                  StepSizeScalingFactor<IV>&) const {}
};
struct SRP_ByConstRef_scale {
  void operator()(StepCount, StepStartAt<IV>,
                  StepSize<IV>&,
                  StepSizeMinAllowedValue<IV>&,
                  const StepSizeScalingFactor<IV>&) const {}
};

// Invalid: rvalue-ref (each param)
struct SRP_ByRvalueRef_dt {
  void operator()(StepCount, StepStartAt<IV>,
                  StepSize<IV>&&,
                  StepSizeMinAllowedValue<IV>&,
                  StepSizeScalingFactor<IV>&) const {}
};
struct SRP_ByRvalueRef_min {
  void operator()(StepCount, StepStartAt<IV>,
                  StepSize<IV>&,
                  StepSizeMinAllowedValue<IV>&&,
                  StepSizeScalingFactor<IV>&) const {}
};
struct SRP_ByRvalueRef_scale {
  void operator()(StepCount, StepStartAt<IV>,
                  StepSize<IV>&,
                  StepSizeMinAllowedValue<IV>&,
                  StepSizeScalingFactor<IV>&&) const {}
};

// Other invalids
struct SRP_ReturnNonVoid {
  int operator()(StepCount, StepStartAt<IV>,
                 StepSize<IV>&, StepSizeMinAllowedValue<IV>&, StepSizeScalingFactor<IV>&) const { return 0; }
};
struct SRP_NonConstOperator {
  void operator()(StepCount, StepStartAt<IV>,
                  StepSize<IV>&, StepSizeMinAllowedValue<IV>&, StepSizeScalingFactor<IV>&) {}
};
struct SRP_WrongIndVar {
  void operator()(StepCount, StepStartAt<float>,
                  StepSize<IV>&, StepSizeMinAllowedValue<IV>&, StepSizeScalingFactor<IV>&) const {}
};
struct SRP_WrongElemTypes_dt {
  void operator()(StepCount, StepStartAt<IV>,
                  StepSize<float>&, StepSizeMinAllowedValue<IV>&, StepSizeScalingFactor<IV>&) const {}
};
struct SRP_WrongElemTypes_min {
  void operator()(StepCount, StepStartAt<IV>,
                  StepSize<IV>&, StepSizeMinAllowedValue<float>&, StepSizeScalingFactor<IV>&) const {}
};
struct SRP_WrongElemTypes_scale {
  void operator()(StepCount, StepStartAt<IV>,
                  StepSize<IV>&, StepSizeMinAllowedValue<IV>&, StepSizeScalingFactor<float>&) const {}
};
struct SRP_Forwarding_dt {
  template<class U>
  void operator()(StepCount, StepStartAt<IV>, U&&,
                  StepSizeMinAllowedValue<IV>&, StepSizeScalingFactor<IV>&) const {}
};
struct SRP_Forwarding_min {
  template<class U>
  void operator()(StepCount, StepStartAt<IV>,
                  StepSize<IV>&, U&&, StepSizeScalingFactor<IV>&) const {}
};
struct SRP_Forwarding_scale {
  template<class U>
  void operator()(StepCount, StepStartAt<IV>,
                  StepSize<IV>&, StepSizeMinAllowedValue<IV>&, U&&) const {}
};
struct SRP_OverloadSet {
  void operator()(StepCount, StepStartAt<IV>,
                  StepSize<IV>&, StepSizeMinAllowedValue<IV>&, StepSizeScalingFactor<IV>&) const {}
  void operator()(StepCount, StepStartAt<IV>,
                  const StepSize<IV>&, StepSizeMinAllowedValue<IV>&, StepSizeScalingFactor<IV>&) const {}
  void operator()(StepCount, StepStartAt<IV>,
                  StepSize<IV>&, const StepSizeMinAllowedValue<IV>&, StepSizeScalingFactor<IV>&) const {}
  void operator()(StepCount, StepStartAt<IV>,
                  StepSize<IV>&, StepSizeMinAllowedValue<IV>&, StepSizeScalingFactor<IV>) const {}
};
struct SRP_NoCall {};

// Single TEST collecting all SRP cases
TEST(StepSizePolicyWithReductionSchemeTrait, AllCases) {
  CHECK_TRAIT2_TRUE(SRP_ByRefVoidConst);
  CHECK_TRAIT2_TRUE(SRP_ByRefVoidConstNoexcept);
  CHECK_TRAIT2_TRUE(SRP_TemplatedExactlyRef);

  CHECK_TRAIT2_FALSE(SRP_ByValue_dt);
  CHECK_TRAIT2_FALSE(SRP_ByValue_min);
  CHECK_TRAIT2_FALSE(SRP_ByValue_scale);

  CHECK_TRAIT2_FALSE(SRP_ByConstRef_dt);
  CHECK_TRAIT2_FALSE(SRP_ByConstRef_min);
  CHECK_TRAIT2_FALSE(SRP_ByConstRef_scale);

  CHECK_TRAIT2_FALSE(SRP_ByRvalueRef_dt);
  CHECK_TRAIT2_FALSE(SRP_ByRvalueRef_min);
  CHECK_TRAIT2_FALSE(SRP_ByRvalueRef_scale);

  CHECK_TRAIT2_FALSE(SRP_ReturnNonVoid);
  CHECK_TRAIT2_FALSE(SRP_NonConstOperator);
  CHECK_TRAIT2_FALSE(SRP_WrongIndVar);
  CHECK_TRAIT2_FALSE(SRP_WrongElemTypes_dt);
  CHECK_TRAIT2_FALSE(SRP_WrongElemTypes_min);
  CHECK_TRAIT2_FALSE(SRP_WrongElemTypes_scale);

  CHECK_TRAIT2_FALSE(SRP_Forwarding_dt);
  CHECK_TRAIT2_FALSE(SRP_Forwarding_min);
  CHECK_TRAIT2_FALSE(SRP_Forwarding_scale);

  CHECK_TRAIT2_FALSE(SRP_OverloadSet);
  CHECK_TRAIT2_FALSE(SRP_NoCall);
}

// ======================= Should be VALID =======================

// Exact signature; other params by value.
struct S_ExactByValue {
  using independent_variable_type = IV;
  using state_type = MyState;
  void operator()(state_type&, StepStartAt<IV>, StepCount, StepSize<IV>) {}
};

// Exact signature; other params by const&.
struct S_ExactConstRefs {
  using independent_variable_type = IV;
  using state_type = MyState;
  void operator()(state_type&,
                  const StepStartAt<IV>&, const StepCount&, const StepSize<IV>&) {}
};

// Same as above but const-qualified and noexcept — should still be valid.
struct S_ConstNoexcept {
  using independent_variable_type = IV;
  using state_type = MyState;
  void operator()(state_type&,
                  const StepStartAt<IV>&, const StepCount&, const StepSize<IV>&) const noexcept {}
};

// ======================= Should be INVALID =======================

// State taken by value → trait forbids being invocable with state&&, which by-value allows.
struct B_StateByValue {
  using independent_variable_type = IV;
  using state_type = MyState;
  void operator()(state_type, StepStartAt<IV>, StepCount, StepSize<IV>) {}
};

// State taken by const& → cannot modify; trait forbids being invocable with const state&.
struct B_StateConstRef {
  using independent_variable_type = IV;
  using state_type = MyState;
  void operator()(const state_type&, StepStartAt<IV>, StepCount, StepSize<IV>) {}
};

// State taken by rvalue ref → cannot bind lvalue state in canonical call.
struct B_StateRvalueRef {
  using independent_variable_type = IV;
  using state_type = MyState;
  void operator()(state_type&&, StepStartAt<IV>, StepCount, StepSize<IV>) {}
};

// Any other param as non-const lvalue ref → canonical call passes prvalues, so not invocable.
struct B_StartAtNonConstLRef {
  using independent_variable_type = IV;
  using state_type = MyState;
  void operator()(state_type&, StepStartAt<IV>&, StepCount, StepSize<IV>) {}
};
struct B_StepCountNonConstLRef {
  using independent_variable_type = IV;
  using state_type = MyState;
  void operator()(state_type&, StepStartAt<IV>, StepCount&, StepSize<IV>) {}
};
struct B_StepSizeNonConstLRef {
  using independent_variable_type = IV;
  using state_type = MyState;
  void operator()(state_type&, StepStartAt<IV>, StepCount, StepSize<IV>&) {}
};

// Any other param as rvalue ref → canonical call OK, but positive check uses const&, which can't bind to &&.
struct B_OtherParamsRvalueRefs {
  using independent_variable_type = IV;
  using state_type = MyState;
  void operator()(state_type&, StepStartAt<IV>&&, StepCount&&, StepSize<IV>&&) {}
};

// Wrong return type.
struct B_ReturnNonVoid {
  using independent_variable_type = IV;
  using state_type = MyState;
  int operator()(state_type&, StepStartAt<IV>, StepCount, StepSize<IV>) { return 0; }
};

// Wrong IV in StepStartAt.
struct B_WrongIndVar {
  using independent_variable_type = IV;
  using state_type = MyState;
  void operator()(state_type&, StepStartAt<float>, StepCount, StepSize<IV>) {}
};

// Wrong StepSize element type.
struct B_WrongStepSizeElem {
  using independent_variable_type = IV;
  using state_type = MyState;
  void operator()(state_type&, StepStartAt<IV>, StepCount, StepSize<float>) {}
};

// Wrong arity / parameter order.
struct B_WrongArity3 {
  using independent_variable_type = IV;
  using state_type = MyState;
  void operator()(state_type&, StepStartAt<IV>, StepCount) {}
};
struct B_WrongOrder {
  using independent_variable_type = IV;
  using state_type = MyState;
  void operator()(state_type&, StepCount, StepStartAt<IV>, StepSize<IV>) {}
};

// Overload set that ALSO accepts const state& → violates "forbid const state&".
struct B_OverloadAlsoConstState {
  using independent_variable_type = IV;
  using state_type = MyState;
  void operator()(state_type&, StepStartAt<IV>, StepCount, StepSize<IV>) {}
  void operator()(const state_type&, StepStartAt<IV>, StepCount, StepSize<IV>) {}
};

// Overload set that ALSO accepts state by value → makes it invocable with state&& (forbidden).
struct B_OverloadAlsoByValue {
  using independent_variable_type = IV;
  using state_type = MyState;
  void operator()(state_type&, StepStartAt<IV>, StepCount, StepSize<IV>) {}
  void operator()(state_type,  StepStartAt<IV>, StepCount, StepSize<IV>) {}
};

// Fully generic forwarding — becomes invocable with const state& and state&& → forbidden.
struct B_TemplatedForwarding {
  using independent_variable_type = IV;
  using state_type = MyState;
  template<class S, class A, class K, class D>
  void operator()(S&&, A&&, K&&, D&&) {}
};

// No operator() at all.
struct B_NoCall {
  using independent_variable_type = IV;
  using state_type = MyState;
};

// =============================== TEST ===============================

#define EXPECT_STEPPER_TRUE(T_)  do { \
  static_assert(PRESSIO_VALUE_OF(StepperWithoutSolver<T_>), #T_ " should be valid"); \
  EXPECT_TRUE((PRESSIO_VALUE_OF(StepperWithoutSolver<T_>))) << #T_; \
} while(0)

#define EXPECT_STEPPER_FALSE(T_) do { \
  static_assert(!PRESSIO_VALUE_OF(StepperWithoutSolver<T_>), #T_ " should be invalid"); \
  EXPECT_FALSE((PRESSIO_VALUE_OF(StepperWithoutSolver<T_>))) << #T_; \
} while(0)

TEST(StepperWithoutSolverTrait, AllCases) {
  // Valid
  EXPECT_STEPPER_TRUE(S_ExactByValue);
  EXPECT_STEPPER_TRUE(S_ExactConstRefs);
  EXPECT_STEPPER_TRUE(S_ConstNoexcept);

  // Invalid
  EXPECT_STEPPER_FALSE(B_StateByValue);
  EXPECT_STEPPER_FALSE(B_StateConstRef);
  EXPECT_STEPPER_FALSE(B_StateRvalueRef);

  EXPECT_STEPPER_FALSE(B_StartAtNonConstLRef);
  EXPECT_STEPPER_FALSE(B_StepCountNonConstLRef);
  EXPECT_STEPPER_FALSE(B_StepSizeNonConstLRef);

  EXPECT_STEPPER_FALSE(B_OtherParamsRvalueRefs);

  EXPECT_STEPPER_FALSE(B_ReturnNonVoid);
  EXPECT_STEPPER_FALSE(B_WrongIndVar);
  EXPECT_STEPPER_FALSE(B_WrongStepSizeElem);
  EXPECT_STEPPER_FALSE(B_WrongArity3);
  EXPECT_STEPPER_FALSE(B_WrongOrder);

  EXPECT_STEPPER_FALSE(B_OverloadAlsoConstState);
  EXPECT_STEPPER_FALSE(B_OverloadAlsoByValue);
  EXPECT_STEPPER_FALSE(B_TemplatedForwarding);
  EXPECT_STEPPER_FALSE(B_NoCall);
}
