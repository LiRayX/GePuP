#pragma once

struct ClassicalRK4
{
  static constexpr int nStages = 4;
  static constexpr double gamma = 0.0;
  static constexpr double c[] = {
      0.0, 0.5, 0.5, 1.0, 1.0};
  static constexpr double aE[][nStages] = {
      {0, 0, 0, 0},
      {0.5, 0, 0, 0},
      {0, 0.5, 0, 0},
      {0, 0, 1, 0},
      {1.0 / 6, 1.0 / 3, 1.0 / 3, 1.0 / 6}};
  static constexpr double aI[][nStages] = {};
};

struct ARK436L2SA
{
  static constexpr int nStages = 6;
  static constexpr double gamma = 1.0 / 4;
  static constexpr double c[] = {
      0.0, 1.0 / 2, 83.0 / 250, 31.0 / 50, 17.0 / 20, 1.0, 1.0};
  static constexpr double aE[][nStages] = {
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.5, 0.0, 0.0, 0.0, 0.0, 0.0},
      {13861.0 / 62500, 6889.0 / 62500, 0.0, 0.0, 0.0, 0.0},
      {-116923316275.0 / 2393684061468, -2731218467317.0 / 15368042101831, 9408046702089.0 / 11113171139209, 0.0, 0.0, 0.0},
      {-451086348788.0 / 2902428689909, -2682348792572.0 / 7519795681897, 12662868775082.0 / 11960479115383, 3355817975965.0 / 11060851509271, 0.0, 0.0},
      {647845179188.0 / 3216320057751, 73281519250.0 / 8382639484533, 552539513391.0 / 3454668386233, 3354512671639.0 / 8306763924573, 4040.0 / 17871, 0.0},
      {82889.0 / 524892, 0.0, 15625.0 / 83664, 69875.0 / 102672, -2260.0 / 8211, 1.0 / 4}};
  static constexpr double aI[][nStages] = {
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {1.0 / 4, 1.0 / 4, 0.0, 0.0, 0.0, 0.0},
      {8611.0 / 62500, -1743.0 / 31250, 1.0 / 4, 0.0, 0.0, 0.0},
      {5012029.0 / 34652500, -654441.0 / 2922500, 174375.0 / 388108, 1.0 / 4, 0.0, 0.0},
      {15267082809.0 / 155376265600, -71443401.0 / 120774400, 730878875.0 / 902184768, 2285395.0 / 8070912, 1.0 / 4, 0.0},
      {82889.0 / 524892, 0.0, 15625.0 / 83664, 69875.0 / 102672, -2260.0 / 8211, 1.0 / 4},
      {82889.0 / 524892, 0.0, 15625.0 / 83664, 69875.0 / 102672, -2260.0 / 8211, 1.0 / 4}};
};
