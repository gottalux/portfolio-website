<template>
  <v-container fluid style="max-width: 1280px">
    <v-row class="text-left">
      <v-col cols="12">
        <h1 class="secondary--text my-3">My Art</h1>
        <v-carousel
          cycle
          height="720"
          hide-delimiters
          next-icon="arrow_forward"
          prev-icon="arrow_back"
        >
          <v-carousel-item
            v-for="(item, i) in art"
            :key="i"
            :src="item.picture"
          >
          </v-carousel-item>
        </v-carousel>
      </v-col>
    </v-row>
  </v-container>
</template>
<script>
export default {
  name: "ResearchPanel",
  data() {
    return {
      art: [
        {
          picture: require("../assets/dackel.png"),
        },
        {
          picture: require("../assets/aChaoticPresentation.png"),
        },
        {
          picture: require("../assets/myUnbrokenHeart.png"),
        },
        {
          picture: require("../assets/forceOfNature.png"),
        },
        {
          picture: require("../assets/mightyPirate.png"),
        },
      ],
    };
  },
  methods: {
    goToUrl(url) {
      window.open(url, "_blank");
    },
  },
};
</script>
