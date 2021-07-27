import React from 'react'
import styled from 'styled-components'

import { ExternalLink } from '@broad/ui'

import DocumentTitle from './DocumentTitle'
import InfoPage from './InfoPage'
import Link from './Link'
import Searchbox from './Searchbox'
import GnomadHeading from './GnomadHeading'

const HomePage = styled(InfoPage)`
  display: flex;
  flex-direction: column;
  align-items: center;
  max-width: 740px;
  margin-top: 90px;
`

const HeadingContainer = styled.div`
  display: flex;
  flex-direction: column;
  align-items: center;
  width: 100%;
  margin-bottom: 40px;
`

const SubHeading = styled.h2`
  padding-top: 0;
  padding-bottom: 0;
  font-size: 1.2em;
  font-weight: lighter;
  letter-spacing: 2px;
  text-align: center;
`

export default () => (
  <HomePage>
    <DocumentTitle />
    <HeadingContainer>
      {/* <GnomadHeading width="60%" /> */}
      <img src="https://ysm-res.cloudinary.com/image/upload/c_fill,f_auto,q_auto:eco,w_920/v1/websites4/live-prod/ycgh/generations/Generations%20Logo_409774_51677_v3.png" width="50%" height="50%"></img>
      <SubHeading>Yale Generations</SubHeading>
      <Searchbox width="100%" />
      <p>
        Examples - Gene:{' '}
        <Link preserveSelectedDataset={false} to="/gene/BRCA2">
          BRCA2
        </Link>
        , Variant:{' '}
        <Link preserveSelectedDataset={false} to="/variant/13-32340300-GT-G">
          13-32340300-GT-G
        </Link>
      </p>
    </HeadingContainer>

    <p> Place holder text
     </p>
  </HomePage>
)
