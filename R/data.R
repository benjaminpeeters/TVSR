#' Monthly policy rates for 137 economies (1946-02 to 2021-07)
#'
#' A data frame of monthly central-bank policy interest rates compiled from
#' public sources (IMF, BIS, national central banks, Wu-Xia shadow series).
#' Used in Peeters, Girard and Gnabo (p2, international monetary spillovers).
#' First three columns are metadata (\code{Countries}, \code{ISO3},
#' \code{`Shadow rates`}); the remaining 903 columns are monthly rates
#' named \code{YYYY-MM} in percentage points (\code{NA} = unavailable).
#'
#' @format A data frame with 137 rows and 906 columns.
#' @source Compiled by the package author from IMF IFS, BIS policy rates
#'   database, national central bank archives, and Wu-Xia (2016) shadow rate
#'   estimates. See p2 supplement, Section A, for the full source list.
"df"

#' Monthly policy rates for 129 economies (policy-only subset)
#'
#' Subset of [df] restricted to economies whose series are based on direct
#' policy rate observations (excluding the Wu-Xia shadow-rate economies).
#' Same column structure as [df].
#'
#' @format A data frame with 129 rows and 906 columns.
#' @source Subset of [df]; see [df] documentation for source details.
"dfPolicyRate"

#' Wu-Xia shadow rates for 8 zero-lower-bound economies
#'
#' Eight economies (EUR, USA, GBR, JPN, CAN, CHE, SWE, DNK) for which the
#' observed policy rate is complemented by the Wu-Xia (2016) shadow rate
#' during zero-lower-bound episodes. Same column structure as [df].
#'
#' @format A data frame with 8 rows and 906 columns.
#' @source Wu and Xia (2016) shadow-rate updates, extended by the package
#'   author using the published estimation code.
"dfPolicyRateShadow"

#' Spatial weight matrix time series for 29 economies (1980-01 to 2019-04)
#'
#' A length-T list of n-by-n spatial weight matrices encoding bilateral
#' trade intensity among 29 economies, used as the \code{W} input to TVSR
#' estimators in the p2 application. Accompanied by the single-period
#' matrix \code{w} (time-average of \code{W}), fixture vectors
#' \code{Cftvw}, \code{sumTrade}, \code{time}, and the cross-sectional
#' output panel \code{Y}.
#'
#' @format \code{W}: a list of length 472 of 29-by-29 numeric matrices
#'   (row-normalised), one per month from 1980-01 to 2019-04.
#'
#'   \code{w}: a 29-by-29 numeric matrix — time-average of \code{W}.
#' @source Bilateral trade data from IMF DOTS and CEPII BACI, aggregated
#'   monthly. See p2 supplement for construction details.
#' @name W
#' @aliases W w
#' @seealso [Y], [time], [sumTrade], [Cftvw]
"W"

#' @rdname W
"w"

#' Output panel: 29 economies, 472 monthly periods
#'
#' The cross-sectional output used as the \code{Y} input to TVSR estimators
#' in the p2 application; monthly policy-rate changes for the 29 economies
#' spanned by [W].
#'
#' @format A 29-by-472 numeric matrix; rows are ISO3 country codes (see
#'   \code{rownames}), columns are monthly dates (see \code{colnames}).
#' @source Derived from [dfPolicyRate]; see p2 supplement.
"Y"

#' POSIXct timestamps aligned with [W] and [Y]
#'
#' @format A length-472 \code{POSIXct} vector covering 1980-01 to 2019-04.
#' @source See [W].
"time"

#' Monthly total bilateral trade volume for the 29-economy panel
#'
#' @format A length-472 numeric vector.
#' @source See [W].
"sumTrade"

#' Monthly trade-intensity factor for the 29-economy panel
#'
#' Scalar per-period trade intensity used to normalise [W]; see the p2
#' supplement for the exact formula.
#'
#' @format A length-472 numeric vector.
#' @source See [W].
"Cftvw"
